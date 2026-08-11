################################################################################
# Build the ICGC Breast-AdenoCa cohort object from the raw tracks
#
# Bins the genome, computes the usable sequence per bin (assembly gaps and the
# ENCODE blacklist removed), averages the eleven covariate tracks onto those
# bins, standardises them, attaches each mutation's covariate values, and
# multiplies copy number by usable sequence to give the exposure the Poisson
# process integrates over.
#
# Usage:  Rscript R/Preprocess_ICGC_BreastAdenoCA.R [tilewidth]
#
#   tilewidth  bin width in bases. Default 2000, the resolution the de novo and
#              refit applications run at. 10000 rebuilds the coarser grid the
#              replication and stability analyses use.
#
# The output is skipped if it already exists, so this is safe to re-run. At
# 2 kb it takes roughly ten minutes and needs about 8 GB.
#
# NOTE on the 10 kb file. The copy already in data/ was built by the predecessor
# project, BEFORE the merge_with_tumor() fix documented in Preprocess_functions.R,
# so rebuilding it here will not reproduce it byte for byte: mutations in
# assembly gaps were previously kept and given an all-zero covariate row, which
# after standardisation reads as an average bin. Rebuilding drops them instead.
# The existing file is left alone unless it is deleted first.
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
TILEWIDTH <- if (length(args)) as.integer(args[1]) else 2000L
stopifnot(TILEWIDTH > 0)

OUT_FILE <- switch(as.character(TILEWIDTH),
                   "2000" = PATH_ICGC2KB,
                   "10000" = PATH_ICGC10KB,
                   file.path(DATA_DIR, sprintf(
                     "ICGC_BreastAdenoCA_avg%dkb_Mutations_Covariates_Copies.rds.gzip",
                     TILEWIDTH %/% 1000L)))

if (file.exists(OUT_FILE)) {
  message("already built: ", basename(OUT_FILE),
          "\n  delete it to rebuild.")
  quit(save = "no", status = 0)
}

# Fail by name now, rather than eight minutes into the binning.
check_inputs(PATHS_PREPROCESS)

message("Building the ", TILEWIDTH, " bp cohort object")
t0 <- Sys.time()
data <- build_icgc_dataset(tilewidth = TILEWIDTH, verbose = TRUE)

# The alignment contract the model depends on, checked before anything is saved.
v <- SignaturePPF_validate(data)
message(sprintf("\n%s mutations | %d samples | %d covariates | %s bins",
                format(v$N, big.mark = ","), v$J, v$p,
                format(v$nbins, big.mark = ",")))
message("covariates: ", paste(colnames(v$SignalTrack), collapse = ", "))

saveRDS(data, OUT_FILE, compress = "gzip")
message(sprintf("\nwrote %s (%.0f MB) in %s", basename(OUT_FILE),
                file.size(OUT_FILE) / 1024^2,
                format(round(Sys.time() - t0))))
