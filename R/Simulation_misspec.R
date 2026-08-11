################################################################################
# Simulation study: behaviour under misspecification
#
# Seven scenarios, each adding one violation of the model's assumptions on top of
# the previous one:
#
#   S0  correctly specified                     (baseline reference)
#   S1  noisy patient-specific epigenome, v^2 = 0.25
#   S2                                  v^2 = 1
#   S3                                  v^2 = 4
#   S4  S3 + hypermutation hotspots for SBS2/SBS13 in 20% of patients
#   S5  S4 + noisy copy-number observations
#   S6  S5 + channel-specific mutation opportunities
#
# 20 replicate datasets per scenario. The LAST 4000 of the 20000 tiles are held
# out of every fit, so every metric is reported in-sample and out-of-sample.
#
# Two use cases are scored:
#   DE NOVO  signatures estimated (K = 15): CompressiveNMF, SignatureAnalyzer,
#            SignaturePPF MAP and MCMC          -> reconstruction metrics
#   FIXED    signatures held at a K = 13 catalogue (the 8 true ones + 5
#            distractors): CompressiveNMF, SignaturePPF MAP and MCMC
#                                               -> attribution and calibration
#
# ---------------------------------------------------------------- RESTARTABILITY
# Every dataset and every fit is skipped if its file is already on disk, so an
# interrupted run resumes where it stopped and a stage can be re-run alone:
#
#   Rscript R/Simulation_misspec.R generate
#   Rscript R/Simulation_misspec.R fit
#   Rscript R/Simulation_misspec.R score
#   Rscript R/Simulation_misspec.R all      # the default
#
# A second argument overrides the number of workers (default 20).
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(foreach)
  library(doParallel)
  library(parallel)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
source(file.path(R_DIR, "Simulation_functions.R"))
source(file.path(R_DIR, "Simulation_functions_misspec.R"))

## ------------------------------------------------------------------ settings
args <- commandArgs(trailingOnly = TRUE)
STAGE <- if (length(args) >= 1) args[1] else "all"
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 20L
stopifnot(STAGE %in% c("all", "generate", "fit", "score"))

J <- 100                 # patients per dataset
N_DATASETS <- 20         # replicates per scenario
N_TEST_BINS <- 4000      # held-out tiles (of 20000)
K_DENOVO <- 15           # upper bound on the number of signatures, de novo
RUN_MCMC <- TRUE

OUT_DIR <- DIR_SIM_MISSPEC

## --------------------------------------------------------------- scenarios
# Cumulative: each scenario adds one violation on top of the previous one.
misspec_scenarios <- function(opportunity = NULL) {
  epi <- 2        # v^2 = 4 epigenome noise, carried from S3 onward
  hot <- list(hotspot_frac = 0.2, hotspot_n = 5, hotspot_mu = 30,
              hotspot_signatures = c("SBS2", "SBS13"))
  cn <- 0.5       # copy-number observation noise, carried from S5 onward
  list(
    "S0_baseline"       = misspec_config(),
    "S1_epigenome_v025" = misspec_config(epigenome_noise_sd = 0.5),
    "S2_epigenome_v1"   = misspec_config(epigenome_noise_sd = 1),
    "S3_epigenome_v4"   = misspec_config(epigenome_noise_sd = 2),
    "S4_hotspots"       = do.call(misspec_config, c(list(epigenome_noise_sd = epi), hot)),
    "S5_cn_noise"       = do.call(misspec_config, c(list(epigenome_noise_sd = epi,
                                                         cn_noise_sd = cn), hot)),
    "S6_opportunity"    = do.call(misspec_config, c(list(epigenome_noise_sd = epi,
                                                         cn_noise_sd = cn,
                                                         opportunity = opportunity), hot))
  )
}

# Computed once, in the PARENT, before any forking: it scans hg19 and caches to
# disk, and twenty workers racing to write the same file would be a bad idea.
opps <- compute_mutation_opportunities_hg19()
scenarios <- misspec_scenarios(opportunity = opps)

## ------------------------------------------------------------------- the grid
# One job per (scenario, replicate). Flattened, so the workers stay busy across
# scenario boundaries rather than draining at the end of each one.
jobs <- expand.grid(replicate = seq_len(N_DATASETS),
                    scenario = names(scenarios),
                    stringsAsFactors = FALSE)
jobs$dir <- file.path(OUT_DIR, jobs$scenario,
                      sprintf("Simulation_%02d", jobs$replicate))
# Deterministic per job, and independent of how the jobs are scheduled.
jobs$seed <- 10000L + 1000L * match(jobs$scenario, names(scenarios)) + jobs$replicate

message(sprintf("%d jobs (%d scenarios x %d replicates) on %d cores | stage: %s",
                nrow(jobs), length(scenarios), N_DATASETS, N_CORES, STAGE))

#' Run one function over the whole job grid, reporting failures by name
#'
#' mclapply returns a try-error rather than throwing, so a job that dies takes
#' the rest of the study with it only if nobody looks. This looks.
run_jobs <- function(label, fun) {
  t0 <- Sys.time()
  out <- mclapply(seq_len(nrow(jobs)), function(i) {
    try(fun(jobs[i, ]), silent = TRUE)
  }, mc.cores = N_CORES, mc.preschedule = FALSE)

  failed <- which(vapply(out, inherits, logical(1), "try-error"))
  if (length(failed)) {
    message(sprintf("\n!! %d of %d %s jobs failed:", length(failed), nrow(jobs), label))
    for (i in failed) {
      message("   ", jobs$scenario[i], "/", sprintf("Simulation_%02d", jobs$replicate[i]),
              ": ", conditionMessage(attr(out[[i]], "condition")))
    }
  }
  message(sprintf("%s: %d/%d ok in %s", label, nrow(jobs) - length(failed), nrow(jobs),
                  format(round(Sys.time() - t0))))
  out
}

################################################################################
# Stage 1 - generate the datasets
################################################################################
if (STAGE %in% c("all", "generate")) {
  message("\n== 1. generate ==")
  invisible(run_jobs("generate", function(job) {
    file <- file.path(job$dir, "data.rds.gzip")
    if (file.exists(file)) return(invisible(NULL))
    create_directory(job$dir)
    set.seed(job$seed, kind = "L'Ecuyer-CMRG")
    # ncores = 1: the parallelism is over replicates, not over patients.
    data <- generate_MutationData_misspec(J = J, ncores = 1,
                                          misspec = scenarios[[job$scenario]])
    saveRDS(data, file, compress = "gzip")
    invisible(NULL)
  }))
}

################################################################################
# Stage 2 - fit the models
################################################################################
if (STAGE %in% c("all", "fit")) {
  message("\n== 2. fit ==")
  invisible(run_jobs("fit", function(job) {
    set.seed(job$seed, kind = "L'Ecuyer-CMRG")
    run_models_misspec(job$dir, K = K_DENOVO, n_test_bins = N_TEST_BINS,
                       run_CompNMF = TRUE, run_SignatureAnalyzer = TRUE,
                       run_MAP = TRUE, run_MCMC = RUN_MCMC,
                       seed = job$seed)
    run_models_misspec_fixed(job$dir, n_test_bins = N_TEST_BINS,
                             run_CompNMF = TRUE,
                             run_MAP = TRUE, run_MCMC = RUN_MCMC,
                             seed = job$seed)
    invisible(NULL)
  }))
}

################################################################################
# Stage 3 - score every fit
################################################################################
if (STAGE %in% c("all", "score")) {
  message("\n== 3. score ==")

  gather <- function(label, fun) {
    out <- run_jobs(label, fun)
    ok <- out[!vapply(out, inherits, logical(1), "try-error")]
    do.call(rbind, ok)
  }

  all_results <- gather("score de novo", function(job)
    postProcessOutput_misspec(job$dir, n_test_bins = N_TEST_BINS))

  fixed_files <- c("output_map_Fixed.rds.gzip", "output_CompNMF_Fixed.rds.gzip",
                   if (RUN_MCMC) "output_mcmc_Fixed.rds.gzip")
  all_calibration <- do.call(rbind, lapply(fixed_files, function(f)
    gather(paste("score", f), function(job)
      postProcessCalibration_misspec(job$dir, fit_file = f,
                                     n_test_bins = N_TEST_BINS))))

  saveRDS(all_results, file.path(OUT_DIR, "all_results.rds"))
  saveRDS(all_calibration, file.path(OUT_DIR, "all_calibration.rds"))
  write.csv(all_results, file.path(OUT_DIR, "all_results.csv"), row.names = FALSE)
  write.csv(all_calibration, file.path(OUT_DIR, "all_calibration.csv"), row.names = FALSE)

  message("\nde novo:  ", nrow(all_results), " rows")
  print(table(all_results$Scenario, all_results$model))
  message("\ncalibration:  ", nrow(all_calibration), " rows")
  print(table(all_calibration$Scenario, all_calibration$model))
}

message("\ndone: outputs in ", OUT_DIR)
