################################################################################
# Build the 80-cancer breast cohort object from the raw calls
#
# Davies et al. (2017). The front end differs from the ICGC one - mutations
# arrive as 80 per-sample CaVEMan VCFs and copy number as 80 ASCAT segment
# tables, so the SNVs have to be read, filtered to clean single-base
# substitutions, and assigned a trinucleotide channel from hg19. Everything
# after that is the same binning code the ICGC preprocessing uses.
#
# Usage:  Rscript R/Preprocess_Breast80.R [tilewidth]
#
#   tilewidth  bin width in bases. Default 10000, which is what the replication
#              analysis needs: it compares this cohort against ICGC on one grid,
#              so the two must be binned identically.
#
# The output is skipped if it already exists. About five minutes.
#
# NOTE. The copy of Breast80_data.rds.gzip shipped in data/ was built by the
# predecessor project, before the merge_with_tumor() fix documented in
# Preprocess_functions.R, so rebuilding will not reproduce it byte for byte:
# mutations in assembly gaps were previously kept and given an all-zero
# covariate row, which after standardisation reads as an average bin. The
# existing file is left alone unless it is deleted first.
#
# (The predecessor's loader hardcoded 10 kb while naming its variables
# `gr_SignalTrack_2kb`. The bins were 10 kb; only the names were wrong.)
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(rtracklayer)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
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
message(sprintf("\n%s mutations | %d samples | %d covariates | %s bins",
                format(v$N, big.mark = ","), v$J, v$p,
                format(v$nbins, big.mark = ",")))
message("covariates: ", paste(colnames(v$SignalTrack), collapse = ", "))

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
