################################################################################
# TensorSignatures comparison, step 1 of 3: build the shared dataset.
#
# Both methods are fitted to the SAME data on the SAME genomic partition - the
# ChromHMM 15-state annotation of breast epithelium (Roadmap E028). The bins ARE
# the ChromHMM segments, so the state each mutation falls in is exact for both
# methods and neither is handicapped by a discretisation the other did not see.
#
# SignaturePPF is given the states as one-hot covariates with `Quies` dropped as
# reference, which makes each beta the log enrichment relative to Quies - the
# same quantity TensorSignatures reports as a state amplitude.
#
# Steps:
#   1. R : this script    - build the dataset, fit SignaturePPF, export the tensor
#   2. sh: bash/run_ts_sweep.sh - fit TensorSignatures over a range of ranks
#   3. R : 03_tensorsignatures_compare.R - import and compare
#
# Runtime: the ChromHMM segmentation is ~600k segments and copy number is built
# per sample, so step 1 takes tens of minutes and a few GB. Everything is cached
# to disk, so rerunning is cheap.
#
# Usage:  Rscript R/02_tensorsignatures_prepare.R
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()

REFERENCE_STATE <- "Quies"
TAG <- "icgc_chromatin"

PATH_DATASET <- file.path(DIR_TENSORSIG, "dataset_chromatin.rds.gzip")
PATH_PPF_FIT <- file.path(DIR_TENSORSIG, "fit_ppf_chromatin.rds.gzip")

################################################################################
# 1. Build the chromatin-state dataset
################################################################################
if (file.exists(PATH_DATASET)) {
  message("using cached dataset: ", basename(PATH_DATASET))
  dat <- readRDS(PATH_DATASET)
} else {
  gr_tumor <- readRDS(PATH_ICGC_SNV)

  # Blacklisted mutations are dropped here rather than relying on bin weights:
  # a mutation inside a blacklisted region has no usable exposure behind it.
  blacklist <- rtracklayer::import(PATH_BLACKLIST)
  gr_tumor <- gr_tumor[-S4Vectors::queryHits(
    GenomicRanges::findOverlaps(gr_tumor, blacklist))]
  gr_tumor <- gr_tumor[GenomicRanges::seqnames(gr_tumor) != "chrY"]

  df_copy <- readr::read_tsv(PATH_ICGC_CN, show_col_types = FALSE)
  gr_copy <- GenomicRanges::GRanges(
    seqnames = paste0("chr", df_copy$chr),
    ranges = IRanges::IRanges(start = df_copy$start, end = df_copy$end),
    strand = "*", sample = df_copy$sampleID, score = df_copy$value)

  dat <- build_chromatin_dataset(gr_tumor, gr_copy, reference = REFERENCE_STATE)
  saveRDS(dat, PATH_DATASET, compress = "gzip")
}

# The dataset must satisfy the model's alignment contract before anything else.
invisible(SignaturePPF_validate(dat))
message("bins: ", nrow(dat$SignalTrack), "  states: ", ncol(dat$SignalTrack) + 1,
        "  samples: ", ncol(dat$CopyTrack))

################################################################################
# 2. Fit SignaturePPF on the chromatin states
#
#    De novo, so that the signature set is estimated rather than assumed - the
#    comparison against TensorSignatures is partly about which signatures each
#    method finds. K is an upper bound; the compressive prior parks the rest.
################################################################################
if (file.exists(PATH_PPF_FIT)) {
  message("using cached PPF fit: ", basename(PATH_PPF_FIT))
  fit <- readRDS(PATH_PPF_FIT)
} else {
  fit <- SignaturePPF(dat,
                      K = 20,
                      method = "map",
                      controls = SignaturePPF_control(maxiter = 500, tol = 1e-6),
                      seed = SEED,
                      verbose = TRUE)
  saveRDS(fit, PATH_PPF_FIT, compress = "gzip")
}
print(fit)

# What the estimated signatures correspond to.
ref_match <- match_to_cosmic(fit$Signatures)
ref_match$mu <- as.numeric(fit$Mu[ref_match$signature])
ref_match <- ref_match[order(-ref_match$mu), ]
write.csv(ref_match, file.path(DIR_TENSORSIG, "ppf_signature_cosmic_match.csv"),
          row.names = FALSE)
print(ref_match)

ggsave(file.path(FIG_DIR, "02_ppf_chromatin_betas.pdf"),
       plot_chromatin_betas(fit, reference = REFERENCE_STATE),
       width = 11, height = 8)

################################################################################
# 3. Export the same data as a TensorSignatures tensor
################################################################################
ts_dir <- export_ts_chromatin(dat, out_dir = DIR_TENSORSIG, tag = TAG)

message("\nNext: fit TensorSignatures with\n",
        "  bash/setup_tensorsig_env.sh      # once\n",
        "  bash/run_ts_sweep.sh             # rank sweep\n",
        "then run R/03_tensorsignatures_compare.R")
