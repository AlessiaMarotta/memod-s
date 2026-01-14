suppressWarnings(
  suppressMessages({
    library(KEGGREST)
    library(fgsea)
    library(data.table)
    library(dplyr)
    library(ggplot2)
  })
)

args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Argomenti mancanti. Uso: Rscript gsea.R <input_dir> <annotation_file> <output_file>")
}

my_dir <- args[1]  # La directory di input
annotation_file <- args[2]
output_file <- args[3]

motif_dirs <- list.dirs(my_dir, recursive = FALSE)

gsea_dir <- file.path(dirname(my_dir), "gsea")
if (!dir.exists(gsea_dir)) {
  dir.create(gsea_dir, recursive = TRUE)
}

emapper.annotations <- read.table(annotation_file, sep = "\t", header = FALSE, quote = "", fill = TRUE, comment.char = "")
if (ncol(emapper.annotations) >= 21) {
    emapper.annotations <- emapper.annotations[, 1:21]
    colnames(emapper.annotations) <- c("query", "seed_ort", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
                                     "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
                                     "CAZy", "BiGG_Reaction", "PFAMs")
} else {
    warning("The annotation file has fewer columns than expected. Please check the format..")
}


get_bed_files <- function(motif_dir) {
  motif_name <- basename(motif_dir)
  subdir <- file.path(motif_dir, "msa")

  if (!dir.exists(subdir)) {
     message(paste("Directory msa non trovata per:", motif_name))
     return(NULL)
  }

  bed_files <- list.files(subdir, pattern = "\\.bed$", full.names = TRUE)

  if (length(bed_files) < 4) {
    warning(paste("Skip motif:", motif_name, "- Only found", length(bed_files), "file BED."))
    return(NULL)
  }

  return(list(
    CDS = as.data.frame(read.table(bed_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    nCDS = as.data.frame(read.table(bed_files[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    tIG = as.data.frame(read.table(bed_files[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    US = as.data.frame(read.table(bed_files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
  ))
}

get_kegg_pathway_info <- function(kegg_code) {
    result <- tryCatch({
      pathway_info <- keggGet(kegg_code)
      if (is.null(pathway_info) || length(pathway_info) == 0) {
        return(list(Name = NA, Description = NA))
      }
      return(list(Name = pathway_info[[1]]$NAME, Description = pathway_info[[1]]$DESCRIPTION))
    }, error = function(e) {
      message("Errore nel recupero di ", kegg_code, ": ", e$message)
      return(list(Name = NA, Description = NA))
    })
    return(result)
}

for (motif_dir in motif_dirs) {
  motif_name <- basename(motif_dir)

  bed_files <- get_bed_files(motif_dir)

  if (is.null(bed_files)) {
    next
  }

  CDS <- bed_files$CDS
  nCDS <- bed_files$nCDS
  tIG <- bed_files$tIG
  US <- bed_files$US

  names(CDS)[names(CDS) == "score"] <- "score_CDS"
  names(nCDS)[names(nCDS) == "score"] <- "score_nCDS"
  names(tIG)[names(tIG) == "score"] <- "score_tIG"
  names(US)[names(US) == "score"] <- "score_US"

  total <- Reduce(function(x, y) merge(x, y, by = c("attribute", "seqid", "start", "end", "gene.annotation"), all = TRUE), list(CDS, nCDS, tIG, US))
  total[is.na(total)] <- 0

  if (nrow(total) == 0) {
      message(paste("Nessun dato unito per", motif_name, "- salto."))
      next
  }

  total$score_tot <- rowSums(total[, c("score_CDS", "score_nCDS", "score_tIG", "score_US")])
  total <- total[order(-total$score_CDS), ]

  total <- total %>%
    mutate(length = abs(end - start))

  if (nrow(total) < 10) {
      message(paste("Troppi pochi punti per LOESS in", motif_name, "- salto."))
      next
  }

  regr <- tryCatch({
      loess(score_CDS ~ length, data = total)
  }, error = function(e) return(NULL))

  if (is.null(regr)) next

  total <- total %>%
    mutate(r_exp = regr$fitted,
           oe = score_CDS - r_exp,
           obs_exp = score_CDS / r_exp)

  genes <- emapper.annotations$query
  KEGG_Pathway <- emapper.annotations$KEGG_Pathway

  total_unique <- total[!duplicated(total$attribute), ]
  ranked_v <- setNames(total_unique$obs_exp, total_unique$attribute)

  ranked_v <- ranked_v[is.finite(ranked_v)]

  pathways <- list()
  for (i in 1:length(KEGG_Pathway)) {
    if (length(KEGG_Pathway) < i || is.na(KEGG_Pathway[i]) || KEGG_Pathway[i] == "") next

    pathways_i <- unlist(strsplit(as.character(KEGG_Pathway[i]), ","))
    for (pathway in pathways_i) {
        pathways[[pathway]] <- c(pathways[[pathway]], genes[i])
    }
  }

  pathways_ko <- pathways[grep("^ko", names(pathways))]

  if (length(pathways_ko) == 0 || length(ranked_v) == 0) {
      message(paste("Dati insufficienti per FGSEA in", motif_name))
      next
  }

  fgseares <- tryCatch({
      fgsea(pathways = pathways_ko,
            stats = ranked_v,
            scoreType = 'std',
            nproc = 1,
            minSize = 1, # Aggiunto per evitare errori su pathway piccoli
            maxSize = 500)
  }, error = function(e) {
      message(paste("Errore FGSEA per", motif_name, ":", e$message))
      return(NULL)
  })

  if (is.null(fgseares) || nrow(fgseares) == 0) next

  combined_pathways <- fgseares %>%
    filter(NES > 0) %>%
    arrange(padj) %>%
    head(16) %>%
    mutate(direction = "Positive") %>%
    bind_rows(
      fgseares %>%
        filter(NES < 0) %>%
        arrange(padj) %>%
        head(16) %>%
        mutate(direction = "Negative")
    )

  if (nrow(combined_pathways) == 0) next

  kegg_codes <- combined_pathways$pathway
  kegg_pathway_info <- sapply(kegg_codes, get_kegg_pathway_info, simplify = FALSE)

  combined_pathways$Name <- sapply(combined_pathways$pathway, function(x) {
    if(x %in% names(kegg_pathway_info)) {
      return(kegg_pathway_info[[x]]$Name)
    } else {
      return(NA)
    }
  })

  combined_pathways <- combined_pathways %>% filter(!is.na(Name))

  if (nrow(combined_pathways) == 0) next

  plot_gsea <- ggplot(combined_pathways, aes(x = reorder(Name, NES), y = NES, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Positive" = "skyblue", "Negative" = "deeppink")) +
    coord_flip() +
    labs(title = paste0("Top significant pathways - ", motif_name), x = "Pathway", y = "Normalized Enrichment Score (NES)") +
    theme_minimal()

  ggsave(filename = file.path(gsea_dir, paste0(motif_name, "_gsea_plot.pdf")), plot = plot_gsea, width = 10, height = 8, dpi = 300)
}

file.create(output_file)
