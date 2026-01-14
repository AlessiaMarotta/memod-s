library(circlize)

args <- commandArgs(TRUE)
my_dir <- args[1]
output_file <- args[2]

motif_dirs <- list.dirs(my_dir, recursive = FALSE)

plot_dir <- file.path(dirname(my_dir), "plot")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

get_bed_files <- function(motif_dir) {
  motif_name <- basename(motif_dir)
  subdir <- file.path(motif_dir, "msa")

  if (!dir.exists(subdir)) {
    warning(paste("msa directory doesn't exist in:", motif_dir))
    return(NULL)
  }

  bed_files <- list.files(subdir, pattern = "\\.bed$", full.names = TRUE)

  cat(paste("Processing:", motif_name, "- Found", length(bed_files), "bed files.\n"))

  if (length(bed_files) < 4) {
    warning(paste("WARNING: < 4 BED files in", subdir, "- Skip this plot."))
    return(NULL)
  }

  return(list(
    my_bed = as.data.frame(read.table(bed_files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    my_bed2 = as.data.frame(read.table(bed_files[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    my_bed3 = as.data.frame(read.table(bed_files[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")),
    my_bed4 = as.data.frame(read.table(bed_files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
  ))
}

for (motif_dir in motif_dirs) {
  motif_name <- basename(motif_dir)

  bed_data <- get_bed_files(motif_dir)

  if (is.null(bed_data)) {
    next
  }

  output_file_path <- file.path(plot_dir, paste0(motif_name, ".pdf"))

  pdf(output_file_path, width = 8, height = 8)

  tryCatch({
      circos.initializeWithIdeogram(bed_data$my_bed4)

      circos.genomicDensity(bed_data$my_bed, col = c("deepskyblue"), track.height = 0.1)
      circos.genomicDensity(bed_data$my_bed2, col = c("deeppink"), track.height = 0.1)
      circos.genomicDensity(bed_data$my_bed3, col = c("forestgreen"), track.height = 0.1)
      circos.genomicDensity(bed_data$my_bed4, col = c("blue"), track.height = 0.1)

      legend("bottomright", legend = c("CDS", "nCDS", "tIG", "US"),
             col = c("deepskyblue", "deeppink", "forestgreen", "blue"),
             pch = 15, pt.cex = 2, cex = 0.8, bty = "n")
  }, error = function(e) {
      message(paste("Error during plotting", motif_name, ":", e$message))
  })

  dev.off()
  circos.clear()
}

file.create(output_file)
