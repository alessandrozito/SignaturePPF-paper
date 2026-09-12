################################################################################
# Produces: no figure. Fits the scenarios behind Figure S8 and Table S5
#
# Sensitivity of the de novo solution to Kmax and to the shrinkage priors
#
#
#   scenario     Kmax   c0    d0     what it varies
#   -----------------------------------------------------------------------
#   K15           15    100   1      Kmax raised
#   K20           20    100   1      Kmax raised further
#   K12_c0_10     12     10   0.1    weaker beta shrinkage
#   K12_c0_1      12      2   0.01   much weaker beta shrinkage
#   K20_c0_10     20     10   0.1    weaker shrinkage, ceiling raised
#   K20_c0_1      20      2   0.01   much weaker shrinkage, ceiling raised
#
#
#
# Usage:  Rscript R/Application_denovo_sensitivity.R [scenario ...]
#
#   scenario   one or more of K15, K20, K12_c0_10, K12_c0_1, K20_c0_10,
#              K20_c0_1. Default is all of them. Given a subset, only those are
#              fitted, which is how two scenarios are put on two cores without
#              them racing for the same output file.
#
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nohup setsid \
#     taskset -c 22 Rscript R/Application_denovo_sensitivity.R \
#     > output/Application_denovo_sensitivity/run.log 2>&1 < /dev/null &
#
# Two scenarios on two cores, which is what the Kmax = 20 pair was run as:
#
#   for s in K20_c0_10 K20_c0_1; do
#     OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nohup setsid \
#       taskset -c $((20 + i)) Rscript R/Application_denovo_sensitivity.R $s \
#       > output/Application_denovo_sensitivity/run_$s.log 2>&1 < /dev/null &
#   done
#
# Each scenario is saved as soon as it finishes and skipped if its file is
# already there, so the run is restartable and a scenario can be added without
# redoing the others.
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()
check_inputs(PATH_ICGC2KB)

MAXITER <- 5000
TOL <- 1e-7

## Which scenarios this process is responsible for. No argument means all of
## them, which is what a single sequential run does.
ALL_SCENARIOS <- c("K15", "K20", "K12_c0_10", "K12_c0_1",
                   "K20_c0_10", "K20_c0_1")
requested <- commandArgs(trailingOnly = TRUE)
if (!length(requested)) requested <- ALL_SCENARIOS
unknown <- setdiff(requested, ALL_SCENARIOS)
if (length(unknown)) {
  stop("unknown scenario(s): ", paste(unknown, collapse = ", "),
       "\n  choose from: ", paste(ALL_SCENARIOS, collapse = ", "), call. = FALSE)
}
wanted <- function(name) name %in% requested

## SEED + 1 is Application_denovo.R's first start. For the two scenarios that
## leave Kmax at 12 that makes the starting point not merely comparable but
## identical to one of the reference's, so any difference is the prior alone.
SEED_START <- SEED + 1L

################################################################################
# Data
################################################################################
data <- readRDS(PATH_ICGC2KB)

message(sprintf("%s mutations | %d samples | %d covariates | %s bins",
                format(length(data$gr_Mutations), big.mark = ","),
                ncol(data$CopyTrack), ncol(data$SignalTrack),
                format(nrow(data$SignalTrack), big.mark = ",")))
message(length(requested), " scenario(s) requested (",
        paste(requested, collapse = ", "), "), one start each, fitted ",
        "sequentially | seed ", SEED_START)

## Printed after each fit, so the log says something useful while the next one
## runs: how many signatures survived the compressive prior, and how long it took.
report <- function(name, fit) {
  message(sprintf("\n>>> %s done: %d of %d signatures above epsilon | %d iterations | %.0f min",
                  name,
                  sum(fit$Mu > 10 * fit$prior$epsilon), ncol(fit$Signatures),
                  as.integer(fit$MAPsolution$iter),
                  as.numeric(fit$runtime, units = "mins")))
}

################################################################################
# K = 15
################################################################################
file_K15 <- file.path(DIR_SENSITIVITY, "MAP_K15.rds.gzip")

if (!wanted("K15")) {
  message("\nK15: not requested, skipping")
} else if (file.exists(file_K15)) {
  message("\nK15: already on disk, skipping")
} else {
  message("\n######## K15 | Kmax = 15, c0 = 100, d0 = 1 ########")
  fit_K15 <- SignaturePPF(
    prune_solution = FALSE,   # the surviving count IS the result; compare later
    data,
    method   = "map",
    sigs     = NULL,
    K        = 15,
    prior    = SignaturePPF_prior(c0 = 100, d0 = 1),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    seed     = SEED_START,
    verbose  = TRUE)
  saveRDS(fit_K15, file_K15, compress = "gzip")
  report("K15", fit_K15)
}

################################################################################
# K = 20
################################################################################
file_K20 <- file.path(DIR_SENSITIVITY, "MAP_K20.rds.gzip")

if (!wanted("K20")) {
  message("\nK20: not requested, skipping")
} else if (file.exists(file_K20)) {
  message("\nK20: already on disk, skipping")
} else {
  message("\n######## K20 | Kmax = 20, c0 = 100, d0 = 1 ########")
  fit_K20 <- SignaturePPF(
    prune_solution = FALSE,   # the surviving count IS the result; compare later
    data,
    method   = "map",
    sigs     = NULL,
    K        = 20,
    prior    = SignaturePPF_prior(c0 = 100, d0 = 1),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    seed     = SEED_START,
    verbose  = TRUE)
  saveRDS(fit_K20, file_K20, compress = "gzip")
  report("K20", fit_K20)
}

################################################################################
# c0 = 10, d0 = 0.1
################################################################################
file_c0_10 <- file.path(DIR_SENSITIVITY, "MAP_K12_c0_10.rds.gzip")

if (!wanted("K12_c0_10")) {
  message("\nK12_c0_10: not requested, skipping")
} else if (file.exists(file_c0_10)) {
  message("\nK12_c0_10: already on disk, skipping")
} else {
  message("\n######## K12_c0_10 | Kmax = 12, c0 = 10, d0 = 0.1 ########")
  fit_c0_10 <- SignaturePPF(
    prune_solution = FALSE,   # the surviving count IS the result; compare later
    data,
    method   = "map",
    sigs     = NULL,
    K        = K_DENOVO,
    prior    = SignaturePPF_prior(c0 = 10, d0 = 0.1),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    seed     = SEED_START,
    verbose  = TRUE)
  saveRDS(fit_c0_10, file_c0_10, compress = "gzip")
  report("K12_c0_10", fit_c0_10)
}

################################################################################
# c0 = 2, d0 = 0.02
################################################################################
file_c0_1 <- file.path(DIR_SENSITIVITY, "MAP_K12_c0_1.rds.gzip")

if (!wanted("K12_c0_1")) {
  message("\nK12_c0_1: not requested, skipping")
} else if (file.exists(file_c0_1)) {
  message("\nK12_c0_1: already on disk, skipping")
} else {
  message("\n######## K12_c0_1 | Kmax = 12, c0 = 2, d0 = 0.01 ########")
  fit_c0_1 <- SignaturePPF(
    prune_solution = FALSE,   # the surviving count IS the result; compare later
    data,
    method   = "map",
    sigs     = NULL,
    K        = K_DENOVO,
    prior    = SignaturePPF_prior(c0 = 2, d0 = 0.01),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    seed     = SEED_START,
    verbose  = TRUE)
  saveRDS(fit_c0_1, file_c0_1, compress = "gzip")
  report("K12_c0_1", fit_c0_1)
}

################################################################################
# Kmax = 20, c0 = 10, d0 = 0.1
################################################################################
file_K20_c0_10 <- file.path(DIR_SENSITIVITY, "MAP_K20_c0_10.rds.gzip")

if (!wanted("K20_c0_10")) {
  message("\nK20_c0_10: not requested, skipping")
} else if (file.exists(file_K20_c0_10)) {
  message("\nK20_c0_10: already on disk, skipping")
} else {
  message("\n######## K20_c0_10 | Kmax = 20, c0 = 10, d0 = 0.1 ########")
  fit_K20_c0_10 <- SignaturePPF(
    prune_solution = FALSE,   # the surviving count IS the result; compare later
    data,
    method   = "map",
    sigs     = NULL,
    K        = 20,
    prior    = SignaturePPF_prior(c0 = 10, d0 = 0.1),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    seed     = SEED_START,
    verbose  = TRUE)
  saveRDS(fit_K20_c0_10, file_K20_c0_10, compress = "gzip")
  report("K20_c0_10", fit_K20_c0_10)
}

################################################################################
# Kmax = 20, c0 = 2, d0 = 0.01
################################################################################
file_K20_c0_1 <- file.path(DIR_SENSITIVITY, "MAP_K20_c0_1.rds.gzip")

if (!wanted("K20_c0_1")) {
  message("\nK20_c0_1: not requested, skipping")
} else if (file.exists(file_K20_c0_1)) {
  message("\nK20_c0_1: already on disk, skipping")
} else {
  message("\n######## K20_c0_1 | Kmax = 20, c0 = 2, d0 = 0.01 ########")
  fit_K20_c0_1 <- SignaturePPF(
    prune_solution = FALSE,   # the surviving count IS the result; compare later
    data,
    method   = "map",
    sigs     = NULL,
    K        = 20,
    prior    = SignaturePPF_prior(c0 = 2, d0 = 0.01),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    seed     = SEED_START,
    verbose  = TRUE)
  saveRDS(fit_K20_c0_1, file_K20_c0_1, compress = "gzip")
  report("K20_c0_1", fit_K20_c0_1)
}

message("\ndone (", paste(requested, collapse = ", "), ") in ", DIR_SENSITIVITY,
        "\n  the comparison against the reference is step 8 of ",
        "R/Reproduce_figures_Application_denovo.R.")
