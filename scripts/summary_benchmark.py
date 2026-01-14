#!/usr/bin/env python3
import pandas as pd
import glob
import os
import sys

bench_dir = snakemake.params.bench_dir
output_raw = snakemake.output.raw
output_stats = snakemake.output.stats

search_pattern = os.path.join(bench_dir, "**", "*_benchmark.txt")
files = glob.glob(search_pattern, recursive=True)

files = [f for f in files if "global_benchmark" not in f]

if not files:
    pd.DataFrame().to_csv(output_raw, sep='\t')
    pd.DataFrame().to_csv(output_stats, sep='\t')
    sys.exit(0)

data_frames = []

for f in files:
    try:
        df = pd.read_csv(f, sep='\t')
        if df.empty: continue

        filename = os.path.basename(f)

        task_name_raw = filename.replace("_benchmark.txt", "")

        sample_name = os.path.basename(os.path.dirname(f))

        if task_name_raw.startswith(sample_name + "_"):
            clean_task = task_name_raw.replace(sample_name + "_", "")
        else:
            clean_task = task_name_raw

        df['Sample'] = sample_name
        df['Task'] = clean_task

        if 'max_rss' in df.columns:
            df['RAM_GB'] = round(df['max_rss'] / 1024, 2)

        if 's' in df.columns:
            df['Time_min'] = round(df['s'] / 60, 2)

        cols = ['Sample', 'Task', 'Time_min', 'RAM_GB', 'mean_load', 'max_rss', 's']
        available = [c for c in cols if c in df.columns]
        data_frames.append(df[available])

    except Exception as e:
        print(f"Error reading {f}: {e}")

if data_frames:
    final_df = pd.concat(data_frames)

    if 'Sample' in final_df.columns:
        final_df = final_df.sort_values(by=['Sample', 'Task'])

    final_df.to_csv(output_raw, sep='\t', index=False)
    print(f"Raw data saved in: {output_raw}")

    if 'Task' in final_df.columns:
        summary_df = final_df.groupby('Task').agg({
            'Time_min': ['mean', 'max', 'std'],
            'RAM_GB': ['mean', 'max']
        }).reset_index()

        summary_df.columns = ['_'.join(col).strip('_') for col in summary_df.columns.values]

        summary_df = summary_df.round(2)

        if 'Time_min_mean' in summary_df.columns:
            summary_df = summary_df.sort_values(by='Time_min_mean', ascending=False)

        summary_df.to_csv(output_stats, sep='\t', index=False)

else:
    open(output_raw, 'w').close()
    open(output_stats, 'w').close()
