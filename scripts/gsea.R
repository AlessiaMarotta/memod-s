library(KEGGREST)
library(fgsea)
library(data.table)
library(dplyr)
library(ggplot2)

args <- commandArgs(TRUE)
my_dir <- args[1]  # La directory di input (dove sono contenute le sottodirectories)
output_file <- args[2]  # File vuoto che verrà creato alla fine
annotation_file <- args[3]

# Trova tutte le sottocartelle all'interno della directory di input
motif_dirs <- list.dirs(my_dir, recursive = FALSE)

# Crea la directory 'gsea' se non esiste
gsea_dir <- file.path(my_dir, "../gsea")
if (!dir.exists(gsea_dir)) {
  dir.create(gsea_dir)
}

# Funzione per ottenere i file BED necessari in ciascuna sottocartella (motivo)
get_bed_files <- function(motif_dir) {
  motif_name <- basename(motif_dir)
  subdir <- file.path(motif_dir, "msa")
  bed_files <- list.files(subdir, pattern = "\\.bed$", full.names = TRUE)

  return(list(
    CDS = as.data.frame(read.table(bed_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    nCDS = as.data.frame(read.table(bed_files[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    tIG = as.data.frame(read.table(bed_files[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    US = as.data.frame(read.table(bed_files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
  ))
}

# Esegui per ogni sottocartella (motivo) trovata
for (motif_dir in motif_dirs) {
  motif_name <- basename(motif_dir)
  bed_files <- get_bed_files(motif_dir)  # Ottieni i file BED
  output_file_path <- file.path(gsea_dir, paste0(motif_name, ".tsv"))}

emapper.annotations <- read.table(args[3], sep = "\t", header = FALSE, quote = "")
colnames(emapper.annotations) <- c("query", "seed_ort", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")

# Assegna i file BED a variabili
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
total$score_tot <- rowSums(total[, c("score_CDS", "score_nCDS", "score_tIG", "score_US")])
total <- total[order(-total$score_tot), ]

total <- total %>%
  mutate(length = abs(end - start))

regr <- loess(score_CDS ~ length, data = total)

total <- total %>%
  mutate(r_exp = regr$fitted,
         oe = score_CDS - r_exp,
         obs_exp = score_CDS / r_exp)

genes <- emapper.annotations$query

KEGG_Pathway <- emapper.annotations$KEGG_Pathway

ranked_v <- setNames(total$obs_exp, total$attribute)

pathways <- list()

for (i in 1:length(KEGG_Pathway)) {
  for (pathway in KEGG_Pathway[[i]]) {
    pathways[[pathway]] <- c(pathways[[pathway]], genes[i])
  }
}

pathways_ko <- pathways[grep("^ko", names(pathways))]

fgseares <- fgsea(pathways = pathways_ko,  # list from all_gene_annotation
                          stats = ranked_v, # ranked list from L
                          scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                          nproc = 1) # for parallelisation
file.create(output_file)

combined_pathways <- fgseares %>%
  filter(pathway != "-") %>%
  filter(NES > 0) %>%
  arrange(padj) %>%
  head(16) %>%
  mutate(direction = "Positive") %>%
  bind_rows(
    fgseares %>%
      filter(pathway != "-") %>%
      filter(NES < 0) %>%
      arrange(padj) %>%
      head(16) %>%
      mutate(direction = "Negative")
  )

kegg_codes <- combined_pathways$pathway

get_kegg_pathway_info <- function(kegg_code) {
  pathway_info <- keggGet(kegg_code)
  if (length(pathway_info) > 0) {
    return(list(Name = pathway_info[[1]]$NAME, Description = pathway_info[[1]]$DESCRIPTION))
  } else {
    return(NA)
  }
}

kegg_pathway_info <- sapply(kegg_codes, get_kegg_pathway_info, simplify = FALSE)

combined_pathways$Name <- sapply(combined_pathways$pathway, function(x) {
  if(x %in% names(kegg_pathway_info)) {
    return(kegg_pathway_info[[x]]$Name)
  } else {
    return(NA)
  }
})

plot_gsea <- ggplot(combined_pathways, aes(x = reorder(Name, NES), y = NES, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Positive" = "skyblue", "Negative" = "deeppink")) +
  coord_flip() +  
  labs(title = "Top significant pathways", x = "Pathway", y = "-log10(padj)") +
  theme_minimal()



# Creazione del plot
#plot_gsea <- ggplot(combined_pathways, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
#  geom_bar(stat = "identity", position = "dodge") +
#  scale_fill_manual(values = c("Positive" = "skyblue", "Negative" = "deeppink")) +
#  coord_flip() +
#  labs(title = paste("Top significant pathways", motif_name), x = "Pathway", y = "-log10(padj)") +
#  theme_minimal()

# Salvataggio del plot
ggsave(filename = file.path(gsea_dir, paste0(motif_name, "_gsea_plot.png")), plot = plot_gsea, width = 8, height = 6, dpi = 300)
