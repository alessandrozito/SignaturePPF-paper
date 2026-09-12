################################################################################
# Produces: no figure. Builds data/preprocessed/ICGC_BreastAdenoCA_avg*.rds.gzip
#
# Build the ICGC Breast-AdenoCa cohort object from the raw tracks
#
# Usage:  Rscript R/Preprocess_ICGC_BreastAdenoCA.R [tilewidth]
#
#   tilewidth  bin width in bases. Default 2000, the resolution the de novo and
#              refit applications run at. 10000 rebuilds the coarser grid the
#              replication and stability analyses use.
#
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
summary_line <- sprintf("\n%s mutations | %d samples | %d covariates | %s bins",
                        format(v$N, big.mark = ","), v$J, v$p,
                        format(v$nbins, big.mark = ","))
covariate_line <- paste(colnames(v$SignalTrack), collapse = ", ")
## v is a processed copy - at 2 kb its reordered CopyTrack alone is 1.2 GB. The
## contract is checked, the numbers are out of it, so let it go before the save.
rm(v); invisible(gc())

message(summary_line)
message("covariates: ", covariate_line)

saveRDS(data, OUT_FILE, compress = "gzip")
message(sprintf("\nwrote %s (%.0f MB) in %s", basename(OUT_FILE),
                file.size(OUT_FILE) / 1024^2,
                format(round(Sys.time() - t0))))
