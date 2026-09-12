################################################################################
# Produces: Figure S1
#
# Main simulation study (Section 4)
#
# Two scenarios differing only in the covariate correlation:
#   Scenario_A_indep   independent covariates
#   Scenario_B_corr    correlated covariates (onion)
#
# 20 replicate datasets each, J = 40 patients, 20000 tiles of 100 bp, 10
# covariates of which the first 5 drive the intensity and the other 5 have
# coefficients of exactly zero.
#
# Eight models. Every one is a variant of PPF except M1 and M7, which are the
# covariate-free competitors:
#
#   M0  MAP, true x        the L = 5 covariates that generated the data, true CN
#   M1  MAP, CompNMF       Poisson NMF, no covariates, no copy number
#   M2  MAP, no x          no covariates, true copy number
#   M3  MAP, all x         all L = 10 covariates, true copy number
#   M4  MCMC, all x        posterior mean for M3, started from its MAP
#   M5  MCMC, all x, D=200 M4 with covariates averaged over 2 consecutive bins
#   M6  MCMC, all x, D=500 M4 with covariates averaged over 5 consecutive bins
#   M7  SignatureAnalyzer  BayesNMF, no covariates, no copy number
#
# Figures:
#   Simualations_results.pdf             the four headline models, labelled (i)-(v)
#   Simualations_results_Supplement2.pdf all eight, labelled M0-M7
#
# Usage:  Rscript R/Simulation_main.R [stage] [cores]
#
#   stage  generate | fit | score | figures | all   (default all)
#   cores  workers, default 20
#
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     taskset -c 0-19 Rscript R/Simulation_main.R
#
# Every dataset and every fit is skipped if already on disk, so an interrupted
# run resumes where it stopped.
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(parallel)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
source(file.path(R_DIR, "Simulation_functions.R"))
source(file.path(R_DIR, "Simulation_functions_misspec.R"))
source(file.path(R_DIR, "Simulation_functions_main.R"))
# For init_from_map(): the chains here start at their own MAP, exactly as the
# applications do, so the same helper serves both.
source(file.path(R_DIR, "Application_functions.R"))

args <- commandArgs(trailingOnly = TRUE)
STAGE <- if (length(args) >= 1) args[1] else "all"
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 20L
stopifnot(STAGE %in% c("all", "generate", "fit", "score", "figures"))

## ------------------------------------------------------------------ settings
J <- 40
N_DATASETS <- 20
K_FIT <- 15                 # upper bound on the number of signatures
DIMS_AGGREG <- c(200, 500)
MAXITER <- 4000
TOL <- 1e-7
NSAMPLES <- 3000
BURNIN <- 1500
SA_K0 <- 15
SA_NRUN <- 5
SA_NITER <- 1e5

SCENARIOS <- c(Scenario_A_indep = "indep", Scenario_B_corr = "onion")
OUT_DIR <- DIR_SIM_MAIN

## ------------------------------------------------------------------- the grid
jobs <- expand.grid(replicate = seq_len(N_DATASETS), scenario = names(SCENARIOS),
                    stringsAsFactors = FALSE)
jobs$dir <- file.path(OUT_DIR, jobs$scenario,
                      sprintf("Simulation_%02d", jobs$replicate))
jobs$seed <- 20000L + 1000L * match(jobs$scenario, names(SCENARIOS)) + jobs$replicate

message(sprintf("%d jobs (%d scenarios x %d replicates) on %d cores | stage: %s",
                nrow(jobs), length(SCENARIOS), N_DATASETS, N_CORES, STAGE))

run_jobs <- function(label, fun) {
  t0 <- Sys.time()
  out <- mclapply(seq_len(nrow(jobs)), function(i) try(fun(jobs[i, ]), silent = TRUE),
                  mc.cores = N_CORES, mc.preschedule = FALSE)
  failed <- which(vapply(out, inherits, logical(1), "try-error"))
  if (length(failed)) {
    message(sprintf("\n!! %d of %d %s jobs failed:", length(failed), nrow(jobs), label))
    for (i in failed) {
      message("   ", jobs$scenario[i], "/",
              sprintf("Simulation_%02d", jobs$replicate[i]), ": ",
              conditionMessage(attr(out[[i]], "condition")))
    }
  }
  message(sprintf("%s: %d/%d ok in %s", label, nrow(jobs) - length(failed),
                  nrow(jobs), format(round(Sys.time() - t0))))
  out
}

################################################################################
# Stage 1 - generate
################################################################################
if (STAGE %in% c("all", "generate")) {
  message("\n== 1. generate ==")
  invisible(run_jobs("generate", function(job) {
    file <- file.path(job$dir, "data.rds.gzip")
    if (file.exists(file)) return(invisible(NULL))
    create_directory(job$dir)
    set.seed(job$seed, kind = "L'Ecuyer-CMRG")
    data <- generate_MutationData(J = J, ncores = 1,
                                  corr = SCENARIOS[[job$scenario]])
    saveRDS(data, file, compress = "gzip")
    invisible(NULL)
  }))
}

################################################################################
# Stage 2 - fit the eight models
################################################################################
fit_all_models <- function(out_dir, seed) {
  data <- readRDS(file.path(out_dir, "data.rds.gzip"))
  tilewidth <- data$simulation_parameters$tilewidth
  Jd <- data$simulation_parameters$J

  SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
  CopyTrackSim <- tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2
  MutMatrix <- SignaturePPF::getTotalMutations(data$gr_Mutations)

  controls_map <- SignaturePPF_control(maxiter = MAXITER, tol = TOL, print_every = 0L)
  controls_mcmc <- SignaturePPF_control(maxiter = MAXITER, tol = TOL,
                                        nsamples = NSAMPLES, burnin = BURNIN,
                                        print_every = 0L)
  prior <- SignaturePPF_prior()

  todo <- function(f) !file.exists(file.path(out_dir, f))
  put <- function(x, f) saveRDS(x, file.path(out_dir, f), compress = "gzip")

  # `prune_solution = FALSE`: the simulation applies its OWN selection rule
  # (Mu above 5 * cutoff AND the spectrum not flat), which is stricter and not
  # the same rule, so the package must hand over every signature it was given
  # and let the scoring decide.
  ppf <- function(dat, ...) SignaturePPF(dat, K = K_FIT, prior = prior,
                                         prune_solution = FALSE,
                                         seed = seed, verbose = FALSE, ...)

  # --- M0: MAP with only the covariates that generated the data ---------------
  if (todo("output_map_TrueCovs.rds.gzip")) {
    covs <- colnames(SignalTrackSim)[seq_len(data$simulation_parameters$p_to_use)]
    gr_true <- data$gr_Mutations
    mcols(gr_true) <- mcols(gr_true)[, c("sample", "channel", covs)]
    put(ppf(list(gr_Mutations = gr_true,
                 SignalTrack = SignalTrackSim[, covs, drop = FALSE],
                 CopyTrack = CopyTrackSim),
            method = "map", controls = controls_map),
        "output_map_TrueCovs.rds.gzip")
  }

  # --- M1: compressive NMF, no covariates, no copy number ---------------------
  if (todo("output_CompNMFBase.rds.gzip")) {
    t0 <- Sys.time()
    out <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = K_FIT, a = prior$a,
                                              alpha = prior$alpha,
                                              epsilon = prior$epsilon, tol = TOL)
    out$time <- Sys.time() - t0
    put(out, "output_CompNMFBase.rds.gzip")
  }

  # --- M2: MAP with copy number but no covariate effect -----------------------
  #     One covariate column has to be present to satisfy the data contract; its
  #     coefficient is pinned at exactly zero, so which one is irrelevant.
  if (todo("output_map_CopyOnly.rds.gzip")) {
    put(ppf(list(gr_Mutations = data$gr_Mutations,
                 SignalTrack = SignalTrackSim, CopyTrack = CopyTrackSim),
            method = "map",
            controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL,
                                            update_Betas = FALSE, print_every = 0L),
            init = SignaturePPF_init(
              Betas_start = matrix(0, ncol(SignalTrackSim), K_FIT))),
        "output_map_CopyOnly.rds.gzip")
  }

  # --- M3: MAP with all covariates -------------------------------------------
  if (todo("output_map_FullModel.rds.gzip")) {
    put(ppf(list(gr_Mutations = data$gr_Mutations,
                 SignalTrack = SignalTrackSim, CopyTrack = CopyTrackSim),
            method = "map", controls = controls_map),
        "output_map_FullModel.rds.gzip")
  }

  # --- M4: MCMC from the M3 mode ---------------------------------------------
  if (todo("output_mcmc_FullModel.rds.gzip")) {
    map <- readRDS(file.path(out_dir, "output_map_FullModel.rds.gzip"))
    put(ppf(list(gr_Mutations = data$gr_Mutations,
                 SignalTrack = SignalTrackSim, CopyTrack = CopyTrackSim),
            method = "mcmc", controls = controls_mcmc,
            init = init_from_map(map, sigs_fixed = FALSE),
            init_mcmc_from_map = FALSE),
        "output_mcmc_FullModel.rds.gzip")
  }

  # --- M5, M6: the same model on coarsened covariates ------------------------
  for (d in DIMS_AGGREG) {
    f_map <- sprintf("output_map_FullModel_agg%d.rds.gzip", d)
    f_mcmc <- sprintf("output_mcmc_FullModel_agg%d.rds.gzip", d)
    if (!todo(f_map) && !todo(f_mcmc)) next
    agg <- aggregate_SignalTrack_CopyTrack(data, lenght_bins = d)
    dat_agg <- list(gr_Mutations = agg$gr_Mutations,
                    SignalTrack = agg$SignalTrack, CopyTrack = agg$CopyTrack)
    if (todo(f_map)) {
      put(ppf(dat_agg, method = "map", controls = controls_map), f_map)
    }
    if (todo(f_mcmc)) {
      map_agg <- readRDS(file.path(out_dir, f_map))
      put(ppf(dat_agg, method = "mcmc", controls = controls_mcmc,
              init = init_from_map(map_agg, sigs_fixed = FALSE),
              init_mcmc_from_map = FALSE), f_mcmc)
    }
  }

  # --- M7: SignatureAnalyzer / BayesNMF --------------------------------------
  if (todo("output_SignatureAnalyzer.rds.gzip")) {
    t0 <- Sys.time()
    out <- sigminer::sig_auto_extract(nmf_matrix = t(MutMatrix), K0 = SA_K0,
                                      nrun = SA_NRUN, niter = SA_NITER,
                                      cores = 1, destdir = tempfile("SA_"))
    out$time <- Sys.time() - t0
    put(out, "output_SignatureAnalyzer.rds.gzip")
  }
  invisible(NULL)
}

if (STAGE %in% c("all", "fit")) {
  message("\n== 2. fit ==")
  invisible(run_jobs("fit", function(job) {
    set.seed(job$seed, kind = "L'Ecuyer-CMRG")
    fit_all_models(job$dir, seed = job$seed)
  }))
}

################################################################################
# Stage 3 - score
################################################################################
# file, model type, and the label the results table carries.
MODEL_FILES <- list(
  list(file = "output_map_TrueCovs.rds.gzip",        type = "PPF",  name = "map_TrueCovs"),
  list(file = "output_CompNMFBase.rds.gzip",         type = "CompNMF", name = "CompNMFBase"),
  list(file = "output_map_CopyOnly.rds.gzip",        type = "PPF",  name = "map_CopyOnly"),
  list(file = "output_map_FullModel.rds.gzip",       type = "PPF",  name = "map_Full"),
  list(file = "output_mcmc_FullModel.rds.gzip",      type = "PPF",  name = "mcmc_Full"),
  list(file = "output_SignatureAnalyzer.rds.gzip",   type = "SignatureAnalyzer",
       name = "SignatureAnalyzer"))

score_one <- function(out_dir) {
  data <- readRDS(file.path(out_dir, "data.rds.gzip"))
  tilewidth <- data$simulation_parameters$tilewidth
  Jd <- data$simulation_parameters$J

  SignalTrackSim <- as.matrix(mcols(data$gr_SignalTrack))[, -1]
  CopyTrackSim <- tilewidth * as.matrix(mcols(data$gr_CopyTrack)) / 2
  CopyTrackNull <- matrix(tilewidth, nrow = nrow(SignalTrackSim), ncol = Jd,
                          dimnames = list(NULL, colnames(CopyTrackSim)))

  Lambda_true <- reconstruct_lambda_raw(SignalTrackSim, CopyTrackSim,
                                        Phi = data$Theta, Betas = data$Betas)
  Theta_total_true <- data$Theta *
    crossprod(exp(SignalTrackSim %*% data$Betas), CopyTrackSim)
  counts_true <- get_counts_from_data(data)

  one <- function(spec) {
    res <- open_rds_file(file.path(out_dir, spec$file))
    if (is.null(res)) return(NULL)
    # The covariate-free competitors have no copy number either, so they are
    # scored against a flat exposure - the only prediction they can make.
    CT <- if (spec$type == "PPF") CopyTrackSim else CopyTrackNull
    est <- getModelEstimates(res, spec$type, SignalTrackSim, CT)
    data.frame(model = spec$name,
               t(compute_RMSE_paramters(res, data, Lambda_true, Theta_total_true,
                                        est$R_hat, est$theta_baseline,
                                        est$Theta_total, est$Beta_hat,
                                        est$Lambda_hat, counts_true)),
               stringsAsFactors = FALSE)
  }
  results <- do.call(rbind, lapply(MODEL_FILES, one))

  # The coarsened models are scored on the FINE grid: their coarse estimates are
  # broadcast back out, so the comparison with the truth is like for like.
  for (d in DIMS_AGGREG) {
    agg <- aggregate_SignalTrack_CopyTrack(data, lenght_bins = d)
    for (m in c("map", "mcmc")) {
      res <- open_rds_file(file.path(out_dir,
                                     sprintf("output_%s_FullModel_agg%d.rds.gzip", m, d)))
      if (is.null(res)) next
      est <- getModelEstimates(res, "PPF", agg$SignalTrackRed, agg$CopyTrackRed)
      results <- rbind(results, data.frame(
        model = sprintf("%s_TrueCovs%d", m, d),
        t(compute_RMSE_paramters(res, data, Lambda_true, Theta_total_true,
                                 est$R_hat, est$theta_baseline, est$Theta_total,
                                 est$Beta_hat, est$Lambda_hat, counts_true)),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(results) && nrow(results)) {
    results$Simulation <- basename(out_dir)
    results$Scenario <- basename(dirname(out_dir))
    results$n <- length(data$gr_Mutations)
  }
  results
}

RESULTS_FILE <- file.path(OUT_DIR, "simulation_results.tsv")

if (STAGE %in% c("all", "score")) {
  message("\n== 3. score ==")
  out <- run_jobs("score", function(job) score_one(job$dir))
  ok <- out[!vapply(out, inherits, logical(1), "try-error")]
  results_all <- do.call(rbind, ok)
  utils::write.table(results_all, RESULTS_FILE, sep = "\t", row.names = FALSE,
                     quote = FALSE)
  message(nrow(results_all), " rows")
  print(table(results_all$Scenario, results_all$model))
}

################################################################################
# Stage 4 - figures
################################################################################
if (STAGE %in% c("all", "figures")) {
  message("\n== 4. figures ==")
  results_all <- utils::read.delim(RESULTS_FILE) |>
    mutate(F1 = 2 * Precision * Sensitivity / (Precision + Sensitivity),
           Scenario2 = ifelse(Scenario == "Scenario_A_indep", "A - indep", "B - corr"))

  # The four headline models, in the order the manuscript lists them.
  MAIN_LABELS <- c(map_TrueCovs      = "(i) MAP, true x",
                   map_Full          = "(ii) MAP, all x",
                   mcmc_Full         = "(iii) MCMC, all x",
                   CompNMFBase       = "(iv) MAP, CompNMF",
                   SignatureAnalyzer = "(v) SignatureAnalyzer")
  MAIN_COLS <- c("#CD2626", "#000D8B", "#89BBF6", "antiquewhite3", "#FF8C00")

  # All eight, labelled as the supplement enumerates them.
  MODEL_LABELS <- c(map_TrueCovs      = "M0. MAP, true x",
                    CompNMFBase       = "M1. MAP, CompNMF",
                    map_CopyOnly      = "M2. MAP, no x",
                    map_Full          = "M3. MAP, all x",
                    mcmc_Full         = "M4. MCMC, all x",
                    mcmc_TrueCovs200  = "M5. MCMC, all x, Delta = 200",
                    mcmc_TrueCovs500  = "M6. MCMC, all x, Delta = 500",
                    SignatureAnalyzer = "M7. SignatureAnalyzer")
  SUPP_COLS <- c("#CD2626", "antiquewhite3", "#FF8C00", "#000D8B", "#89BBF6",
                 "lightgreen", "forestgreen", "#9467BD")

  # Beta is undefined for the models that carry no covariates, so it is left
  # blank rather than drawn as the norm of the true effects.
  NO_BETAS <- c("CompNMFBase", "map_CopyOnly", "SignatureAnalyzer")

  panel_data <- function(labels) {
    results_all |>
      filter(.data$model %in% names(labels)) |>
      mutate(Model = factor(labels[.data$model], levels = unname(labels)),
             Signatures = .data$rmse_sig,
             Theta = .data$rmse_theta,
             Betas = ifelse(.data$model %in% NO_BETAS, NA, .data$rmse_Betas),
             Lambda = .data$rmse_lambda)
  }

  ## ---------------------------------------------------------- main figure
  p_main <- panel_data(MAIN_LABELS) |>
    dplyr::select("Model", "Scenario2", "Signatures", "Theta", "Betas", "Lambda") |>
    tidyr::gather("key", "value", -"Model", -"Scenario2") |>
    mutate(key = factor(.data$key,
                        levels = c("Signatures", "Theta", "Betas", "Lambda"))) |>
    ggplot() +
    geom_boxplot(aes(x = .data$Scenario2, y = log(.data$value),
                     fill = .data$Model, colour = .data$Model), alpha = 0.6) +
    facet_wrap(~ key, scales = "free", nrow = 1) +
    scale_fill_manual(name = "Model", values = MAIN_COLS) +
    scale_colour_manual(name = "Model", values = MAIN_COLS) +
    theme_bw() + theme(aspect.ratio = 1) +
    xlab("Scenario") + ylab("log RMSE")
  ggsave(file.path(FIG_DIR, "Simualations_results.pdf"), p_main,
         width = 10.6, height = 2.34)

  ## ---------------------------------------------------------- supplement
  p_supp <- panel_data(MODEL_LABELS) |>
    mutate(K = .data$Kest,
           Lambda = log(.data$Lambda), Signatures = log(.data$Signatures),
           Theta = log(.data$Theta), Betas = log(.data$Betas)) |>
    dplyr::select("Model", "Scenario2", "K", "F1", "Lambda", "Signatures",
                  "Theta", "Betas") |>
    tidyr::gather("key", "value", -"Model", -"Scenario2") |>
    mutate(key = factor(.data$key, levels = c("K", "F1", "Lambda", "Signatures",
                                              "Theta", "Betas"))) |>
    ggplot() +
    geom_boxplot(aes(x = .data$Scenario2, y = .data$value,
                     fill = .data$Model, colour = .data$Model), alpha = 0.6) +
    facet_wrap(~ key, scales = "free", nrow = 2) +
    scale_fill_manual(name = "Model", values = SUPP_COLS) +
    scale_colour_manual(name = "Model", values = SUPP_COLS) +
    theme_bw() + theme(aspect.ratio = 1, axis.title.y = element_blank()) +
    xlab("Scenario")
  ggsave(file.path(FIG_DIR, "Simualations_results_Supplement2.pdf"), p_supp,
         width = 9.9, height = 5.20)

  ## ------------------------------------------- Tables S1 and S2, as CSV
  summary_tbl <- results_all |>
    filter(.data$model %in% names(MODEL_LABELS)) |>
    mutate(Model = MODEL_LABELS[.data$model]) |>
    group_by(.data$Model, .data$Scenario2) |>
    summarise(across(c("time", "iter", "Kest", "F1", "rmse_sig", "rmse_theta",
                       "rmse_Betas", "rmse_lambda", "effectiveBetas",
                       "effectiveSigs", "effectiveTheta", "effectiveMu",
                       "effectiveSigma2", "effectiveLogPost"),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~stats::sd(.x, na.rm = TRUE)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop")
  write.csv(summary_tbl, file.path(OUT_DIR, "simulation_summary.csv"),
            row.names = FALSE)
  print(as.data.frame(summary_tbl[, c("Model", "Scenario2", "Kest_mean",
                                      "F1_mean", "time_mean", "iter_mean")]))
}

message("\ndone: outputs in ", OUT_DIR)
