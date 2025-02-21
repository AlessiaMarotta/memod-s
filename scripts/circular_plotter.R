library(circlize)

# Ottieni il percorso della directory di input e il file vuoto da creare
args <- commandArgs(TRUE)
my_dir <- args[1]  # La directory di input (dove sono contenuti i sottodirectory/motivi)
output_file <- args[2]  # File vuoto che verrà creato alla fine

# Trova tutte le sottocartelle all'interno della directory di input
motif_dirs <- list.dirs(my_dir, recursive = FALSE)
#motif_dirs <- motif_dirs[grep("msa$", motif_dirs)]

# Crea la directory 'plot' se non esiste

plot_dir <- file.path(dirname(my_dir), "plot")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}


# Funzione per ottenere i file BED necessari in ciascuna sottocartella (motivo)
get_bed_files <- function(motif_dir) {
  motif_name <- basename(motif_dir)
  subdir <- file.path(motif_dir, "msa")
  bed_files <- list.files(subdir, pattern = "\\.bed$", full.names = TRUE)
  
  return(list(
    my_bed = as.data.frame(read.table(bed_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    my_bed2 = as.data.frame(read.table(bed_files[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    my_bed3 = as.data.frame(read.table(bed_files[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    my_bed4 = as.data.frame(read.table(bed_files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
  ))
}

# Esegui il plot per ogni sottocartella (motivo) trovata
for (motif_dir in motif_dirs) {
  motif_name <- basename(motif_dir)
  bed_files <- get_bed_files(motif_dir)  # Ottieni i file BED
  
  # Crea il percorso per il file PDF di output
  output_file_path <- file.path(plot_dir, paste0(motif_name, ".pdf"))
  
  # Inizializzazione del plot circolare
  pdf(output_file_path, width = 8, height = 8)
  
  circos.initializeWithIdeogram(bed_files$my_bed4)
  
  # Specifica della densità di metilazione per ciascun BED
  circos.genomicDensity(bed_files$my_bed, col = c("deepskyblue"), track.height = 0.1)
  circos.genomicDensity(bed_files$my_bed2, col = c("deeppink"), track.height = 0.1)
  circos.genomicDensity(bed_files$my_bed3, col = c("forestgreen"), track.height = 0.1)
  circos.genomicDensity(bed_files$my_bed4, col = c("blue"), track.height = 0.1)
  
  legend("bottomright", legend = c("CDS", "nCDS", "tIG", "US"),
         col = c("deepskyblue", "deeppink", "forestgreen", "blue"),
         pch = 15, pt.cex = 2, cex = 0.8, bty = "n")
  
  dev.off()
}

# Crea un file vuoto con il nome del file specificato come output
file.create(output_file)
