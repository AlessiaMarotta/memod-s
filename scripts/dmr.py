#!/usr/bin/env python3
"""
dmr_sliding_segmentation.py

Sliding windows (step=1 bp) su tutto il genoma + segmentazione + volcano + dotplot.

Input:
 - control TSV: must contain columns Contig, Position, Modified_bases, Unmodified_bases
 - treated TSV: same format
 - contig_lengths TSV: Contig <tab> Length (header optional)

Outputs:
 - window-level results (TSV)
 - segment-level results (TSV)
 - volcano plot (PNG)  -> saved to --output
 - dotplot (PNG)      -> saved to --output with suffix _dotplot.png

Requisiti opzionali:
 - ruptures per la segmentazione (pip install ruptures)
"""
import argparse
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import sys
import os

# try import ruptures (se non presente, il codice segnala e procede senza segmentazione)
try:
    import ruptures as rpt
except Exception:
    rpt = None


def benjamini_hochberg(pvals):
    """Return adjusted p-values (FDR) using BH procedure.
    Input: array-like p-values (floats). Output: numpy array of adj p-values in [0,1].
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return np.array([])
    # Handle NaN by setting them to 1 so they don't affect ordering
    nan_mask = np.isnan(p)
    p_safe = p.copy()
    p_safe[nan_mask] = 1.0
    order = np.argsort(p_safe)
    ranks = np.empty(n, int)
    ranks[order] = np.arange(1, n + 1)
    adj = p_safe * n / ranks
    # enforce monotonicity: adj_corrected[i] = min_{j>=i} adj[j]
    adj_corrected = np.minimum.accumulate(adj[::-1])[::-1]
    adj_corrected[adj_corrected > 1] = 1.0
    # restore NaNs as NaN
    adj_corrected[nan_mask] = np.nan
    return adj_corrected


def aggregate_segment(df_segment, eps=1e-6):
    """Aggrega i counts e calcola statistiche su un dataframe di windows contigui (segmento)."""
    mod_ctrl = int(df_segment["Modified_bases_ctrl"].sum())
    unmod_ctrl = int(df_segment["Unmodified_bases_ctrl"].sum())
    mod_treat = int(df_segment["Modified_bases_treat"].sum())
    unmod_treat = int(df_segment["Unmodified_bases_treat"].sum())
    total_ctrl = mod_ctrl + unmod_ctrl
    total_treat = mod_treat + unmod_treat
    pct_ctrl = mod_ctrl / (total_ctrl + eps)
    pct_treat = mod_treat / (total_treat + eps)
    delta = pct_treat - pct_ctrl
    # fisher
    try:
        _, p = fisher_exact([[mod_treat, unmod_treat], [mod_ctrl, unmod_ctrl]])
    except Exception:
        p = 1.0
    # log2 fold change of percentages (small offset to avoid div0)
    log2fc = np.log2((pct_treat + eps) / (pct_ctrl + eps))
    ratio = (mod_treat - mod_ctrl) / (mod_treat + mod_ctrl + eps)
    return {
        "Modified_bases_ctrl": mod_ctrl,
        "Unmodified_bases_ctrl": unmod_ctrl,
        "Modified_bases_treat": mod_treat,
        "Unmodified_bases_treat": unmod_treat,
        "Percent_ctrl": pct_ctrl,
        "Percent_treat": pct_treat,
        "delta_percent": delta,
        "p_value": p,
        "log2FC": log2fc,
        "ratio_diff": ratio
    }


def parse_contig_lengths(path):
    """Legge contig_lengths TSV/CSV; supporta header oppure file due colonne senza header."""
    try:
        df = pd.read_csv(path, sep="\t", header=0, comment='#', engine='python')
        cols = [c.lower() for c in df.columns]
        if "contig" in cols and ("length" in cols or "len" in cols):
            rename_map = {}
            for c in df.columns:
                if c.lower() == "contig":
                    rename_map[c] = "Contig"
                if c.lower() in ("length", "len"):
                    rename_map[c] = "Length"
            df = df.rename(columns=rename_map)[["Contig", "Length"]]
            return df
    except Exception:
        pass
    # fallback: due colonne senza header
    df = pd.read_csv(path, sep="\t", header=None, comment='#', engine='python')
    if df.shape[1] >= 2:
        df = df.iloc[:, :2]
        df.columns = ["Contig", "Length"]
        return df
    raise ValueError("Impossibile leggere contig_lengths. Fornisci Contig <tab> Length (header accettato).")


def ensure_cols(df):
    """Verifica e normalizza colonne input sample (case-insensitive)."""
    cols_map = {c.lower(): c for c in df.columns}
    needed = ["contig", "position", "modified_bases", "unmodified_bases"]
    for n in needed:
        if n not in cols_map:
            raise ValueError(f"Colonna richiesta '{n}' non trovata nel file. Colonne presenti: {list(df.columns)}")
    df = df.rename(columns={
        cols_map["contig"]: "Contig",
        cols_map["position"]: "Position",
        cols_map["modified_bases"]: "Modified_bases",
        cols_map["unmodified_bases"]: "Unmodified_bases"
    })
    # assicura tipi interi
    df["Position"] = df["Position"].astype(int)
    df["Modified_bases"] = df["Modified_bases"].astype(int)
    df["Unmodified_bases"] = df["Unmodified_bases"].astype(int)
    return df


def sliding_counts_for_contig(contig, length, ctrl_df, treat_df, window):
    """
    Costruisce vettori di sliding-sum per questa contig:
    - ctrl_mod_sliding: somma Modified_bases nelle finestre sliding
    - ctrl_unmod_sliding...
    ritorna dict con arrays e start_bp array
    """
    L = int(length)
    if L < window:
        return None  # niente sliding windows possibili
    # arrays per posizione (0-based indexing: posizione-1 => index)
    ctrl_mod_pos = np.zeros(L, dtype=np.int64)
    ctrl_unmod_pos = np.zeros(L, dtype=np.int64)
    treat_mod_pos = np.zeros(L, dtype=np.int64)
    treat_unmod_pos = np.zeros(L, dtype=np.int64)

    # riempi gli array con i counts per posizione (attenzione a pos > length)
    # assumiamo che Position sia 1-based. convertiamo a 0-based index
    for _, r in ctrl_df.iterrows():
        pos = int(r["Position"]) - 1
        if 0 <= pos < L:
            ctrl_mod_pos[pos] += int(r["Modified_bases"])
            ctrl_unmod_pos[pos] += int(r["Unmodified_bases"])
    for _, r in treat_df.iterrows():
        pos = int(r["Position"]) - 1
        if 0 <= pos < L:
            treat_mod_pos[pos] += int(r["Modified_bases"])
            treat_unmod_pos[pos] += int(r["Unmodified_bases"])

    # sliding sum efficienti: convolution con kernel di 1 di lunghezza window
    kernel = np.ones(window, dtype=np.int64)
    mod_ctrl_sl = np.convolve(ctrl_mod_pos, kernel, mode='valid')
    unmod_ctrl_sl = np.convolve(ctrl_unmod_pos, kernel, mode='valid')
    mod_treat_sl = np.convolve(treat_mod_pos, kernel, mode='valid')
    unmod_treat_sl = np.convolve(treat_unmod_pos, kernel, mode='valid')

    starts = np.arange(0, L - window + 1, dtype=np.int64)
    ends = starts + window - 1

    return {
        "starts": starts,
        "ends": ends,
        "mod_ctrl": mod_ctrl_sl,
        "unmod_ctrl": unmod_ctrl_sl,
        "mod_treat": mod_treat_sl,
        "unmod_treat": unmod_treat_sl
    }


def main():
    parser = argparse.ArgumentParser(description="DMR: sliding windows (step=1) + segmentation + volcano + dotplot")
    parser.add_argument("--control", required=True, help="TSV control (Contig, Position, Modified_bases, Unmodified_bases)")
    parser.add_argument("--treated", required=True, help="TSV treated")
    parser.add_argument("--contig_lengths", required=True, help="TSV with Contig and Length")
    parser.add_argument("--window", type=int, default=50, help="Window size in bp (default 50)")
    parser.add_argument("--output", required=True, help="Output PNG file for the volcano plot (dotplot saved as <this>_dotplot.png)")
    parser.add_argument("--results", required=True, help="Output TSV with window-level results")
    parser.add_argument("--segments_results", default=None, help="Output TSV with segment-level results (default: <results>_segments.tsv)")
    parser.add_argument("--seg_method", choices=["pelt", "binseg"], default="pelt", help="Segmentation method (ruptures): pelt or binseg")
    parser.add_argument("--seg_penalty", type=float, default=5.0, help="Penalty for PELT (higher -> fewer breakpoints). Tune per dataset.")
    parser.add_argument("--min_cov", type=int, default=10, help="Minimum combined coverage (ctrl+treat) in a window to include it in segmentation (default:10)")
    parser.add_argument("--min_windows", type=int, default=3, help="Minimum number of windows per final segment (small segments will be merged) (default:3)")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff to call significant DMR (default: 0.05)")
    parser.add_argument("--eps", type=float, default=1e-6, help="Small epsilon to avoid div0")
    args = parser.parse_args()

    if args.segments_results is None:
        args.segments_results = args.results.rstrip(".tsv") + "_segments.tsv"

    if rpt is None:
        print("ATTENZIONE: la libreria 'ruptures' non è installata. La segmentazione sarà saltata.", file=sys.stderr)

    # Carica dati
    ctrl = pd.read_csv(args.control, sep="\t", comment='#', engine='python')
    treat = pd.read_csv(args.treated, sep="\t", comment='#', engine='python')
    ctrl = ensure_cols(ctrl)
    treat = ensure_cols(treat)
    contig_df = parse_contig_lengths(args.contig_lengths)
    contig_df["Length"] = contig_df["Length"].astype(int)

    # Prepara accumulo risultati windows (vettoriale per contig)
    win_rows = []

    # Processa contig per contig (riduce uso RAM)
    for _, crow in contig_df.iterrows():
        contig = crow["Contig"]
        length = int(crow["Length"])
        print(f"Processing contig {contig}, length {length} bp ...")
        c_ctrl = ctrl[ctrl["Contig"] == contig]
        c_treat = treat[treat["Contig"] == contig]

        sliding = sliding_counts_for_contig(contig, length, c_ctrl, c_treat, args.window)
        if sliding is None:
            print(f"  Contig {contig} shorter than window ({length} < {args.window}). Skipping sliding windows.")
            continue

        starts = sliding["starts"]
        ends = sliding["ends"]
        mod_ctrl_arr = sliding["mod_ctrl"].astype(np.int64)
        unmod_ctrl_arr = sliding["unmod_ctrl"].astype(np.int64)
        mod_treat_arr = sliding["mod_treat"].astype(np.int64)
        unmod_treat_arr = sliding["unmod_treat"].astype(np.int64)

        # calcoli vettoriali: percentuali, p-values, log2fc, ratio
        total_ctrl_arr = mod_ctrl_arr + unmod_ctrl_arr
        total_treat_arr = mod_treat_arr + unmod_treat_arr
        pct_ctrl_arr = mod_ctrl_arr / (total_ctrl_arr + args.eps)
        pct_treat_arr = mod_treat_arr / (total_treat_arr + args.eps)
        delta_pct_arr = pct_treat_arr - pct_ctrl_arr
        pvals = np.ones_like(pct_ctrl_arr, dtype=float)
        log2fc_arr = np.log2((pct_treat_arr + args.eps) / (pct_ctrl_arr + args.eps))
        ratio_arr = (mod_treat_arr - mod_ctrl_arr) / (mod_treat_arr + mod_ctrl_arr + args.eps)
        neglog_p = np.zeros_like(pvals)

        # calcolo p-values finestra-per-finestra
        for i in range(len(pvals)):
            a = int(mod_treat_arr[i])
            b = int(unmod_treat_arr[i])
            c = int(mod_ctrl_arr[i])
            d = int(unmod_ctrl_arr[i])
            try:
                _, p = fisher_exact([[a, b], [c, d]])
            except Exception:
                p = 1.0
            pvals[i] = p
            neglog_p[i] = -np.log10(p + 1e-300)

        # append a win_rows
        for i in range(len(starts)):
            win_rows.append({
                "Contig": contig,
                "Window_index": int(starts[i]),  # start coordinate as index
                "start_bp": int(starts[i]),
                "end_bp": int(ends[i]),
                "Modified_bases_ctrl": int(mod_ctrl_arr[i]),
                "Unmodified_bases_ctrl": int(unmod_ctrl_arr[i]),
                "Modified_bases_treat": int(mod_treat_arr[i]),
                "Unmodified_bases_treat": int(unmod_treat_arr[i]),
                "total_ctrl": int(total_ctrl_arr[i]),
                "total_treat": int(total_treat_arr[i]),
                "tot_cov": int(total_ctrl_arr[i] + total_treat_arr[i]),
                "Percent_ctrl": float(pct_ctrl_arr[i]),
                "Percent_treat": float(pct_treat_arr[i]),
                "delta_percent": float(delta_pct_arr[i]),
                "p_value": float(pvals[i]),
                "log2FC": float(log2fc_arr[i]),
                "ratio_diff": float(ratio_arr[i]),
                "neg_log10_p": float(neglog_p[i])
            })

    # DataFrame windows
    win_df = pd.DataFrame(win_rows)
    if win_df.empty:
        print("Nessuna window generata. Termino.", file=sys.stderr)
        sys.exit(1)

    # ordina e salva
    win_df = win_df.sort_values(["Contig", "start_bp"]).reset_index(drop=True)
    win_df.to_csv(args.results, sep="\t", index=False)
    print(f"Saved window-level results to {args.results} (n_windows = {len(win_df)})")

    # Se ruptures non disponibile, salviamo win_df e terminiamo (ma generiamo anche il dotplot basato su segmenti vuoti)
    if rpt is None:
        print("ruptures non disponibile: salto segmentazione. Genero solo plots placeholder.", file=sys.stderr)
        empty_seg = pd.DataFrame(columns=[
            "Contig", "start_window", "end_window", "start_bp", "end_bp", "n_windows",
            "Modified_bases_ctrl", "Unmodified_bases_ctrl", "Modified_bases_treat", "Unmodified_bases_treat",
            "Percent_ctrl", "Percent_treat", "delta_percent", "p_value", "log2FC"
        ])
        empty_seg.to_csv(args.segments_results, sep="\t", index=False)
        print(f"Saved empty segments file to {args.segments_results}")
        # placeholder dotplot
        dotplot_path = args.output.replace(".png", "_dotplot.png")
        plt.figure(figsize=(6, 3))
        plt.text(0.5, 0.5, "ruptures not installed\nno segmentation performed", ha="center", va="center")
        plt.axis("off")
        plt.savefig(dotplot_path)
        print(f"Saved placeholder dotplot to {dotplot_path}")
        sys.exit(0)

    # ---------- Segmentazione per contig (usando ratio_diff) ----------
    segments = []
    contigs = win_df["Contig"].unique()
    for contig in contigs:
        dfc = win_df[win_df["Contig"] == contig].copy().reset_index(drop=True)
        if dfc.empty:
            continue
        # mask di copertura minima
        mask_cov = dfc["tot_cov"] >= args.min_cov
        n = len(dfc)
        i = 0
        while i < n:
            if not mask_cov.iloc[i]:
                i += 1
                continue
            j = i
            while j < n and mask_cov.iloc[j]:
                j += 1
            region = dfc.iloc[i:j].reset_index(drop=True)
            L = len(region)
            if L == 0:
                i = j
                continue
            if L < args.min_windows:
                agg = aggregate_segment(region, args.eps)
                segments.append({
                    "Contig": contig,
                    "start_window": int(region["Window_index"].min()),
                    "end_window": int(region["Window_index"].max()),
                    "start_bp": int(region["start_bp"].min()),
                    "end_bp": int(region["end_bp"].max()),
                    "n_windows": L,
                    **agg
                })
                i = j
                continue

            signal = region["ratio_diff"].values.reshape(-1, 1)
            if args.seg_method == "pelt":
                algo = rpt.Pelt(model="l2").fit(signal)
                bkps = algo.predict(pen=args.seg_penalty)
            else:
                max_bkps = max(1, int(L / max(10, args.seg_penalty)))
                algo = rpt.Binseg(model="l2").fit(signal)
                bkps = algo.predict(n_bkps=max_bkps)

            prev = 0
            for b in bkps:
                seg_region = region.iloc[prev:b]
                if len(seg_region) == 0:
                    prev = b
                    continue
                agg = aggregate_segment(seg_region, args.eps)
                segments.append({
                    "Contig": contig,
                    "start_window": int(seg_region["Window_index"].min()),
                    "end_window": int(seg_region["Window_index"].max()),
                    "start_bp": int(seg_region["start_bp"].min()),
                    "end_bp": int(seg_region["end_bp"].max()),
                    "n_windows": int(len(seg_region)),
                    **agg
                })
                prev = b
            i = j

    seg_df = pd.DataFrame(segments)
    if seg_df.empty:
        print("Nessun segmento generato (nessuna regione con copertura sufficiente).", file=sys.stderr)
        seg_df.to_csv(args.segments_results, sep="\t", index=False)
        sys.exit(0)

    # Merge e pulizia segmenti troppo piccoli (come nel tuo script originale)
    seg_df = seg_df.sort_values(["Contig", "start_window"]).reset_index(drop=True)

    def merge_two_rows(df, idx_a, idx_b):
        a = df.loc[idx_a]
        b = df.loc[idx_b]
        new = {}
        new["Contig"] = a["Contig"]
        new["start_window"] = min(a["start_window"], b["start_window"])
        new["end_window"] = max(a["end_window"], b["end_window"])
        new["start_bp"] = min(a["start_bp"], b["start_bp"])
        new["end_bp"] = max(a["end_bp"], b["end_bp"])
        new["n_windows"] = int((new["end_window"] - new["start_window"]) + 1)
        for col in ["Modified_bases_ctrl", "Unmodified_bases_ctrl", "Modified_bases_treat", "Unmodified_bases_treat"]:
            new[col] = int(a[col] + b[col])
        total_ctrl = new["Modified_bases_ctrl"] + new["Unmodified_bases_ctrl"]
        total_treat = new["Modified_bases_treat"] + new["Unmodified_bases_treat"]
        new["Percent_ctrl"] = new["Modified_bases_ctrl"] / (total_ctrl + args.eps)
        new["Percent_treat"] = new["Modified_bases_treat"] / (total_treat + args.eps)
        new["delta_percent"] = new["Percent_treat"] - new["Percent_ctrl"]
        try:
            _, p = fisher_exact([[new["Modified_bases_treat"], new["Unmodified_bases_treat"]], [new["Modified_bases_ctrl"], new["Unmodified_bases_ctrl"]]])
        except Exception:
            p = 1.0
        new["p_value"] = p
        new["log2FC"] = np.log2((new["Percent_treat"] + args.eps) / (new["Percent_ctrl"] + args.eps))
        return pd.Series(new)

    # ripeti merging come nel tuo script
    while True:
        small_idx = seg_df[(seg_df["n_windows"] < args.min_windows)].index
        if len(small_idx) == 0 or len(seg_df) == 1:
            break
        idx = small_idx[0]
        if idx == 0:
            neighbor = idx + 1
        elif idx == len(seg_df) - 1:
            neighbor = idx - 1
        else:
            cov_prev = seg_df.loc[idx - 1, "Modified_bases_ctrl"] + seg_df.loc[idx - 1, "Unmodified_bases_ctrl"] + seg_df.loc[idx - 1, "Modified_bases_treat"] + seg_df.loc[idx - 1, "Unmodified_bases_treat"]
            cov_next = seg_df.loc[idx + 1, "Modified_bases_ctrl"] + seg_df.loc[idx + 1, "Unmodified_bases_ctrl"] + seg_df.loc[idx + 1, "Modified_bases_treat"] + seg_df.loc[idx + 1, "Unmodified_bases_treat"]
            neighbor = idx - 1 if cov_prev >= cov_next else idx + 1
        a = min(idx, neighbor)
        b = max(idx, neighbor)
        new_row = merge_two_rows(seg_df, a, b)
        top = seg_df.iloc[:a]
        bottom = seg_df.iloc[b + 1:]
        seg_df = pd.concat([top, pd.DataFrame([new_row]), bottom], ignore_index=True)
        seg_df = seg_df.sort_values(["Contig", "start_window"]).reset_index(drop=True)

    # Ricalcola FDR e metriche (robusto)
    adj = benjamini_hochberg(seg_df["p_value"].values)
    # evitiamo FDR esattamente 0: clamp inferiore (es. 1e-16) per evitare -log10 -> infinito
    min_adj = 1e-16
    # adj può contenere NaN (se p_value era NaN): manteniamo NaN
    adj_clipped = np.where(np.isnan(adj), np.nan, np.clip(adj, min_adj, 1.0))
    seg_df["adj_p"] = adj_clipped
    seg_df["FDR_significant"] = seg_df["adj_p"] < args.fdr
    seg_df["neg_log10_p"] = -np.log10(seg_df["p_value"] + 1e-300)
    # y-axis basata su FDR (adj p-values)
    seg_df["neg_log10_adj_p"] = -np.log10(seg_df["adj_p"] + 1e-300)
    seg_df["abs_log2FC"] = np.abs(seg_df["log2FC"])

    # salva segmenti
    seg_df.to_csv(args.segments_results, sep="\t", index=False)
    print(f"Saved segment-level results to {args.segments_results}")
    print(f"Identified {len(seg_df)} segments across {len(contigs)} contigs. {seg_df['FDR_significant'].sum()} segments FDR-significant at {args.fdr}")

    # ---------- Volcano plot (segmenti) con FDR sull'asse y ----------
    volcano_path = args.output
    plt.figure(figsize=(9, 6))
    ax = plt.gca()
    ycol = "neg_log10_adj_p"

    # punti non significativi
    ax.scatter(seg_df["log2FC"], seg_df[ycol], c="grey", alpha=0.6, label="segments")

    # punti significativi: FDR_significant + abs(log2FC)>1
    sig_mask = seg_df["FDR_significant"].fillna(False) & (seg_df["abs_log2FC"] > 1)
    if sig_mask.any():
        ax.scatter(seg_df.loc[sig_mask, "log2FC"], seg_df.loc[sig_mask, ycol],
                   c="red", alpha=0.9, label="FDR-significant & |log2FC|>1")

    # linea soglia per FDR
    ax.axhline(-np.log10(args.fdr), color="blue", linestyle="--", label=f"FDR = {args.fdr}")

    # linee verticali per cutoff di log2FC
    ax.axvline(1, color="green", linestyle="--")
    ax.axvline(-1, color="green", linestyle="--")

    ax.set_xlabel("log2 Fold Change (Treated vs Control) [segments]")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"Volcano plot (segment-level DMRs) - penalty={args.seg_penalty}, min_cov={args.min_cov}, min_windows={args.min_windows}")
    ax.legend(loc="upper right")

    # --- secondary axis: mostra il valore FDR corrispondente (10^(-y)) sulla destra ---
    #try:
    #    secax = ax.secondary_yaxis('right', functions=(lambda y: 10 ** (-y), lambda f: -np.log10(f)))
    #    secax.set_ylabel('FDR')
    #    import matplotlib.ticker as mticker
    #    secax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1e'))
        # rimuovi eventuale offset text (es. 1e+02) che a volte compare in basso a destra
    #    secax.yaxis.get_offset_text().set_visible(False)
    #except Exception:
        # ambiente matplotlib molto vecchio o mancante di secondary_yaxis: ignoriamo
    #    pass

    plt.tight_layout()
    plt.savefig(volcano_path)
    plt.close()
    print(f"Saved volcano plot to {volcano_path}")

    # ---------- Dot plot segmenti ----------
    epsilon = args.eps  # per evitare divisione per zero
    seg_df["meth_diff_ratio"] = (seg_df["Modified_bases_treat"] - seg_df["Modified_bases_ctrl"]) / \
                                (seg_df["Modified_bases_treat"] + seg_df["Modified_bases_ctrl"] + epsilon)
    # per il colore uso -log10(adj_p) (FDR) * segno della differenza, è più interpretabile
    seg_df["color_metric"] = -np.log10(seg_df["adj_p"] + 1e-300) * np.sign(seg_df["meth_diff_ratio"]) * np.abs(seg_df["meth_diff_ratio"])
    vmax = np.abs(seg_df["color_metric"]).max() if len(seg_df) > 0 else 1.0

    plt.figure(figsize=(14, 4))
    ax = plt.gca()
    import matplotlib as mpl
    cmap = mpl.cm.bwr
    norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
    for _, row in seg_df.iterrows():
        color = cmap(norm(row["color_metric"]))
        ax.hlines(
            y=row["meth_diff_ratio"],
            xmin=row["start_bp"],
            xmax=row["end_bp"],
            colors=[color],
            linewidth=3
        )
    ax.axhline(0, color="black", linestyle="--")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylabel("(Sample2 - Sample1) / (Sample2 + Sample1)")
    ax.set_title("Segment-level methylation difference ratio")
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("-log10(FDR) * |ratio| (signed by sign(ratio))")
    plt.tight_layout()
    dotplot_path = volcano_path.replace(".png", "_dotplot.png")
    plt.savefig(dotplot_path)
    plt.close()
    print(f"Saved dot plot to {dotplot_path}")


if __name__ == "__main__":
    main()
