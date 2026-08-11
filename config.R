## Paths and shared settings for the analysis pipeline.
##
## Sourced by every script. Nothing here has side effects beyond creating the
## output directories, so it is safe to source repeatedly.
##
## Every path derives from one root, so the pipeline moves to a cluster by
## setting one environment variable:
##
##   SIGNATUREPPF_PAPER   this repository (default ~/SignaturePPF-paper)
##
## All inputs live in data/ inside the repository, so the analyses depend on
## nothing outside it. SIGNATUREPPF_DATA overrides that location if the data has
## to sit on a different volume - on a cluster it is usually too large for a home
## directory quota.

PAPER_ROOT <- Sys.getenv("SIGNATUREPPF_PAPER",
                         unset = path.expand("~/SignaturePPF-paper"))
DATA_DIR   <- Sys.getenv("SIGNATUREPPF_DATA",
                         unset = file.path(PAPER_ROOT, "data"))

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
## E028 is the Roadmap ChromHMM segmentation of breast epithelium.
PATH_CHROMHMM  <- file.path(DATA_DIR, "E028_15_coreMarks_dense.bed")
PATH_BLACKLIST <- file.path(DATA_DIR, "hg19-blacklist.v2.bed")
PATH_GAPS      <- file.path(DATA_DIR, "gaps_hg19.bed")
PATH_ICGC_SNV  <- file.path(DATA_DIR, "Breast-AdenoCa_snp.rds.gzip")
PATH_ICGC_CN   <- file.path(DATA_DIR, "20170119_final_consensus_copynumber_donor")

## Fail early and by name, rather than three steps into a pipeline.
check_inputs <- function(paths = c(PATH_BREAST80, PATH_ICGC10KB, PATH_CHROMHMM,
                                   PATH_BLACKLIST, PATH_GAPS, PATH_ICGC_SNV,
                                   PATH_ICGC_CN)) {
  gone <- paths[!file.exists(paths)]
  if (length(gone)) {
    stop("missing input file(s):\n  ", paste(gone, collapse = "\n  "),
         "\n\nSee data/README.md.", call. = FALSE)
  }
  invisible(TRUE)
}

## --------------------------------------------------------------- output dirs
DIR_REPLICATION <- file.path(OUTPUT_DIR, "Replication_80Breast")
DIR_TENSORSIG   <- file.path(OUTPUT_DIR, "Comparison_TensorSignatures")
DIR_STABILITY   <- file.path(OUTPUT_DIR, "Covariate_stability")
DIR_SIM_MISSPEC <- file.path(OUTPUT_DIR, "Simulation_misspec")
for (d in c(DIR_REPLICATION, DIR_TENSORSIG, DIR_STABILITY, DIR_SIM_MISSPEC)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## The 96 hg19 trinucleotide opportunities, computed once by the simulation
## study and cached. A derived quantity, so it lives in output/ rather than
## data/ - nothing outside the pipeline is needed to rebuild it.
PATH_OPPORTUNITY <- file.path(OUTPUT_DIR, "mutation_opportunities_hg19.rds")

## ------------------------------------------------------------------ analysis
## The reference signatures refitted in the replication analysis. Chosen in the
## exploratory pass as the COSMIC signatures with support in breast cohorts, plus
## the MMRd set (SBS6/20/26/44) the 80-cohort was selected to contain.
SIGS_TO_USE <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                 "SBS20", "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30",
                 "SBS40a", "SBS44")

SEED <- 10L

## ------------------------------------------------------------- BLAS threads
BLAS_THREADS <- as.integer(Sys.getenv("SIGNATUREPPF_BLAS_THREADS", unset = "1"))
if (BLAS_THREADS > 0) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(BLAS_THREADS)
  } else {
    message("RhpcBLASctl is not installed, so BLAS threads are left at the ",
            "system default.\n  install.packages(\"RhpcBLASctl\"), or start R ",
            "with OPENBLAS_NUM_THREADS=", BLAS_THREADS, ".")
  }
}

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
