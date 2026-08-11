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

## ------------------------------------------------------- covariate bigWigs
## The eleven genomic covariates, as the tracks they are built from. The seven
## histone/CTCF marks come in a tissue and a cell-line version and are averaged;
## the other four have a single source each.
PATH_GC        <- file.path(DATA_DIR, "gc_content_1kb.bigWig")
PATH_METHYL    <- file.path(DATA_DIR, "Breast-Cancer_Methylation.bigWig")
PATH_REPLITIME <- file.path(DATA_DIR, "wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig")
PATH_NUCLEOSOME <- file.path(DATA_DIR,
                             "GSM920557_hg19_wgEncodeSydhNsomeK562Sig_1kb.bigWig")

CHROMATIN_MARKS <- c("CTCF", "H3K9me3", "H3K36me3", "H3K27me3", "H3K27ac",
                     "H3K4me1", "H3K4me3")
path_mark <- function(mark, source = c("tissue", "cell")) {
  source <- match.arg(source)
  file.path(DATA_DIR, sprintf("Breast-Cancer_%s_%s_2kb.bigWig", source, mark))
}

## The preprocessed cohort at the 2 kb resolution the applications run at, built
## by R/Preprocess_ICGC_BreastAdenoCA.R from everything above.
PATH_ICGC2KB <- file.path(DATA_DIR,
                          "ICGC_BreastAdenoCA_avg2kb_Mutations_Covariates_Copies.rds.gzip")

## The covariate tracks and the two masks - shared by both cohorts.
PATHS_TRACKS <- c(PATH_BLACKLIST, PATH_GAPS,
                  PATH_GC, PATH_METHYL, PATH_REPLITIME, PATH_NUCLEOSOME,
                  vapply(CHROMATIN_MARKS, path_mark, "", source = "tissue"),
                  vapply(CHROMATIN_MARKS, path_mark, "", source = "cell"))

## Everything the ICGC preprocessing reads.
PATHS_PREPROCESS <- c(PATH_ICGC_SNV, PATH_ICGC_CN, PATHS_TRACKS)

## ------------------------------------------------ the 80-cancer cohort, raw
## Davies et al. (2017): per-sample CaVEMan calls and ASCAT copy-number
## segments, one file each per tumour.
PATH_BREAST80_SNV <- file.path(DATA_DIR, "SNP80Breast")
PATH_BREAST80_CN  <- file.path(DATA_DIR, "copyNumber80Breast")
PATHS_PREPROCESS_BREAST80 <- c(PATH_BREAST80_SNV, PATH_BREAST80_CN, PATHS_TRACKS)

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
DIR_SIM_MAIN    <- file.path(OUTPUT_DIR, "Simulation_main")
DIR_DENOVO      <- file.path(OUTPUT_DIR, "Application_denovo")
DIR_REFIT       <- file.path(OUTPUT_DIR, "Application_refit")
DIR_SENSITIVITY <- file.path(OUTPUT_DIR, "Application_denovo_sensitivity")
for (d in c(DIR_REPLICATION, DIR_TENSORSIG, DIR_STABILITY, DIR_SIM_MISSPEC,
            DIR_SIM_MAIN,
            DIR_DENOVO, DIR_REFIT, DIR_SENSITIVITY)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## The 96 hg19 trinucleotide opportunities, computed once by the simulation
## study and cached. A derived quantity, so it lives in output/ rather than
## data/ - nothing outside the pipeline is needed to rebuild it.
PATH_OPPORTUNITY <- file.path(OUTPUT_DIR, "mutation_opportunities_hg19.rds")

## ------------------------------------------------------------------ analysis
## The reference catalogue held fixed wherever signatures are not estimated: the
## main refit application AND the two-cohort replication. Chosen in the
## exploratory pass as the COSMIC signatures with support in breast cohorts, plus
## the MMRd set (SBS6/20/26/44) the 80-cancer cohort was selected to contain.
##
## One list for both on purpose. The compressive prior still applies, so a
## signature the cohort does not support is parked near epsilon rather than
## forced onto the data - which means handing both analyses the same fifteen
## costs nothing and makes their relevance weights directly comparable.
SIGS_TO_USE <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS6", "SBS8",
                 "SBS20", "SBS26", "SBS17a", "SBS17b", "SBS18", "SBS30",
                 "SBS40a", "SBS44")

## Upper bound on the number of signatures in the de novo application.
K_DENOVO <- 12

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
                    "Preprocess_functions.R",
                    "Application_functions.R",
                    "TensorSignatures_functions.R")

load_functions <- function() {
  for (f in FUNCTION_FILES) source(file.path(R_DIR, f))
  invisible(TRUE)
}
