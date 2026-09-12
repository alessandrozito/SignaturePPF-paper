################################################################################
# Produces: no figure. Builds data/preprocessed/Breast80_data.rds.gzip
#
# Build the 80-cancer breast cohort object from the raw calls
#
#
# Usage:  Rscript R/Preprocess_Breast80.R [tilewidth]
#
#   tilewidth  bin width in bases. Default 10000, which is what the replication
#              analysis needs: it compares this cohort against ICGC on one grid,
#              so the two must be binned identically.
#
# The output is skipped if it already exists. About five minutes.
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(rtracklayer)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()

args <- commandArgs(trailingOnly = TRUE)
TILEWIDTH <- if (length(args)) as.integer(args[1]) else 10000L
stopifnot(TILEWIDTH > 0)

OUT_FILE <- if (TILEWIDTH == 10000L) PATH_BREAST80 else
  file.path(DATA_DIR, sprintf("Breast80_data_%dkb.rds.gzip", TILEWIDTH %/% 1000L))

if (file.exists(OUT_FILE)) {
  message("already built: ", basename(OUT_FILE), "\n  delete it to rebuild.")
  quit(save = "no", status = 0)
}

check_inputs(PATHS_PREPROCESS_BREAST80)

message("Building the ", TILEWIDTH, " bp cohort object")
t0 <- Sys.time()
data <- build_breast80_dataset(tilewidth = TILEWIDTH, verbose = TRUE)

v <- SignaturePPF_validate(data)
summary_line <- sprintf("\n%s mutations | %d samples | %d covariates | %s bins",
                        format(v$N, big.mark = ","), v$J, v$p,
                        format(v$nbins, big.mark = ","))
covariate_line <- paste(colnames(v$SignalTrack), collapse = ", ")
## v is a processed copy of the whole cohort. The contract is checked and the
## numbers are out of it, so let it go - the ICGC object gets loaded just below
## for the grid comparison, and there is no reason to hold both.
rm(v); invisible(gc())

message(summary_line)
message("covariates: ", covariate_line)

# The replication analysis puts this cohort and ICGC side by side, so the two
# grids have to agree. Checked here rather than discovered three minutes into
# that script.
if (file.exists(PATH_ICGC10KB) && TILEWIDTH == 10000L) {
  icgc <- readRDS(PATH_ICGC10KB)
  same_bins <- nrow(icgc$SignalTrack) == nrow(data$SignalTrack)
  same_covs <- identical(colnames(icgc$SignalTrack), colnames(data$SignalTrack))
  message(sprintf("\nvs the ICGC 10 kb object: same bin count %s | same covariates %s",
                  same_bins, same_covs))
  if (!same_bins || !same_covs) {
    warning("this cohort and the ICGC 10 kb object are not on the same grid; ",
            "the replication analysis will refuse to run.", call. = FALSE)
  }
}

saveRDS(data, OUT_FILE, compress = "gzip")
message(sprintf("\nwrote %s (%.0f MB) in %s", basename(OUT_FILE),
                file.size(OUT_FILE) / 1024^2,
                format(round(Sys.time() - t0))))
