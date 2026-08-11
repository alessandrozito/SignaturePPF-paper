## Paths and shared settings for the analysis pipeline.
##
## Sourced by every script. Nothing here has side effects beyond creating the
## output directories, so it is safe to source repeatedly.
##
## Every path is derived from one of two roots, so the pipeline can be moved or
## run on a cluster by setting two environment variables and nothing else:
##
##   SIGNATUREPPF_PAPER   this repository            (default ~/SignaturePPF-paper)
##   SIGNATUREPPF_DATA    preprocessed cohort data   (default ~/SigPoisProcess/data)
##
## SIGNATUREPPF_DATA still points into the old project because the preprocessed
## `data` objects were built there and the preprocessing has not been ported to
## SignaturePPF yet. Once SignaturePPF_preprocess() exists, those objects get
## rebuilt into data/ here and the default changes.

PAPER_ROOT <- Sys.getenv("SIGNATUREPPF_PAPER",
                         unset = path.expand("~/SignaturePPF-paper"))
DATA_DIR   <- Sys.getenv("SIGNATUREPPF_DATA",
                         unset = path.expand("~/SigPoisProcess/data"))

OUTPUT_DIR <- file.path(PAPER_ROOT, "output")
FIG_DIR    <- file.path(PAPER_ROOT, "figures")
R_DIR      <- file.path(PAPER_ROOT, "R")

for (d in c(OUTPUT_DIR, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- input data
## Preprocessed cohort objects: lists of gr_Mutations / SignalTrack / CopyTrack.
## Both are on the SAME 10 kb bin grid and carry the SAME 11 covariates, which is
## what makes the replication comparison meaningful.
PATH_BREAST80 <- file.path(DATA_DIR, "Breast80_data.rds.gzip")
PATH_ICGC10KB <- file.path(DATA_DIR,
                           "ICGC_BreastAdenoCA_avg10kb_Mutations_Covariates_Copies.rds.gzip")

## Raw inputs for the chromatin-state comparison.
PATH_CHROMHMM  <- Sys.getenv("SIGNATUREPPF_CHROMHMM",
                             unset = path.expand("~/E028_15_coreMarks_dense.bed"))
PATH_BLACKLIST <- file.path(DATA_DIR, "data_for_application/hg19-blacklist.v2.bed")
PATH_GAPS      <- file.path(DATA_DIR, "data_for_application/gaps_hg19.bed")
PATH_ICGC_SNV  <- file.path(DATA_DIR, "data_for_application/Breast-AdenoCa_snp.rds.gzip")
PATH_ICGC_CN   <- file.path(DATA_DIR,
                            "data_for_application/20170119_final_consensus_copynumber_donor")

## --------------------------------------------------------------- output dirs
DIR_REPLICATION <- file.path(OUTPUT_DIR, "Replication_80Breast")
DIR_TENSORSIG   <- file.path(OUTPUT_DIR, "Comparison_TensorSignatures")
for (d in c(DIR_REPLICATION, DIR_TENSORSIG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## ------------------------------------------------------------------ analysis
## The reference signatures refitted in the replication analysis. Chosen in the
## exploratory pass as the COSMIC signatures with support in breast cohorts, plus
## the MMRd set (SBS6/20/26/44) the 80-cohort was selected to contain.
SIGS_TO_USE <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                 "SBS20", "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30",
                 "SBS40a", "SBS44")

SEED <- 10L

## --------------------------------------------------------- TensorSignatures
## The conda environment built by setup_tensorsig_env.sh. TensorSignatures
## 0.5.0 pins tensorflow <= 1.15, whose wheels stop at cp37, so it cannot share an
## interpreter with anything modern - see that script.
TENSORSIG_PYTHON <- Sys.getenv(
  "TENSORSIG_PYTHON",
  unset = path.expand("~/miniconda3/envs/tensorsig/bin/python"))

## The helper files, listed explicitly: R/ is flat, so sourcing everything in it
## would also source the analysis scripts.
FUNCTION_FILES <- c("Utils_functions.R",
                    "Plot_functions.R",
                    "TensorSignatures_functions.R")

load_functions <- function() {
  for (f in FUNCTION_FILES) source(file.path(R_DIR, f))
  invisible(TRUE)
}
