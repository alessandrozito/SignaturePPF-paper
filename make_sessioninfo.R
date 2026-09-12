## Load everything the pipeline touches, so the recorded versions are complete.
pkgs <- c("SignaturePPF", "tidyverse", "patchwork", "GenomicRanges", "rtracklayer",
          "BSgenome", "BSgenome.Hsapiens.UCSC.hg19", "Biostrings", "GenomicFeatures",
          "ggalluvial", "corrplot", "RcppHungarian", "RhpcBLASctl", "readxl",
          "cluster", "Rcpp")
loaded <- character(); failed <- character()
for (p in pkgs) {
  ok <- suppressWarnings(suppressPackageStartupMessages(
          require(p, character.only = TRUE, quietly = TRUE)))
  if (ok) loaded <- c(loaded, p) else failed <- c(failed, p)
}
source("config.R"); load_functions()

hdr <- c(
  "Session information for the SignaturePPF-paper pipeline",
  paste0("Generated ", format(Sys.Date())),
  "",
  "Every package the analysis scripts use was loaded, then config.R sourced,",
  "before calling sessionInfo(). Regenerate with:",
  "",
  "    Rscript make_sessioninfo.R",
  "")
if (length(failed)) hdr <- c(hdr, paste("NOT INSTALLED:", paste(failed, collapse = ", ")), "")
writeLines(c(hdr, capture.output(sessionInfo())), "SESSIONINFO.txt")
cat("loaded", length(loaded), "of", length(pkgs), "\n")
if (length(failed)) cat("not installed:", paste(failed, collapse = ", "), "\n")
