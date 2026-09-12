################################################################################
# Produces: no figure. Fits the fixed-signature model, which
#           Figures 4 and S6 are drawn from
#
# Refit application: ICGC Breast-AdenoCa at 2 kb, signatures held fixed
#
# The fifteen COSMIC signatures with support in breast cohorts are held at their
# catalogue values and only the activities and the covariate coefficients are
# estimated.
#
#   MAP     one run - with the signatures fixed the objective is far better
#           behaved than the de novo one, so multiple starts buy nothing
#   MCMC    10000 iterations, 5000 burn-in, started at the MAP
#
# Usage:  Rscript R/Application_refit.R
#
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nohup setsid \
#     taskset -c 1 Rscript R/Application_refit.R \
#     > output/Application_refit/run.log 2>&1 < /dev/null &
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(ggplot2)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()
check_inputs(PATH_ICGC2KB)

## ------------------------------------------------------------------ settings
MAXITER <- 5000
TOL <- 1e-7
NSAMPLES <- 10000
BURNIN <- 5000
LOGPOST_EVERY <- 10

################################################################################
# 1. Data and the reference catalogue
################################################################################
data <- readRDS(PATH_ICGC2KB)
message(sprintf("%s mutations | %d samples | %d covariates | %s bins",
                format(length(data$gr_Mutations), big.mark = ","),
                ncol(data$CopyTrack), ncol(data$SignalTrack),
                format(nrow(data$SignalTrack), big.mark = ",")))

CosmicSigs <- COSMIC_v3.4_SBS96_GRCh37[, SIGS_TO_USE]
message("refitting ", ncol(CosmicSigs), " fixed signatures: ",
        paste(SIGS_TO_USE, collapse = ", "))

prior <- SignaturePPF_prior()

################################################################################
# 2. MAP
################################################################################
message("\n== MAP ==")
map_file <- file.path(DIR_REFIT, "MAPSolution.rds.gzip")
if (file.exists(map_file)) {
  message("using existing MAP: ", basename(map_file))
  map <- readRDS(map_file)
} else {
  # NOT pruned: this mode is the chain's starting point and has to stay at the
  # full K (the package rejects a narrower R_start).
  map <- SignaturePPF(data,
                      prune_solution = FALSE,
                      sigs = CosmicSigs, sigs_fixed = TRUE,
                      method = "map",
                      prior = prior,
                      controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
                      seed = SEED, verbose = TRUE)
  saveRDS(map, map_file, compress = "gzip")
}
print(map)

################################################################################
# 3. MCMC, started from the MAP above. Single starting point
################################################################################
message("\n== MCMC ==")
message(sprintf("starting the chain from the MAP, log posterior %s",
                format(map_logposterior(map), big.mark = ",", nsmall = 2)))

mcmc_file <- file.path(DIR_REFIT, "MCMCSolution.rds.gzip")

if (file.exists(mcmc_file)) {
  message("using the existing chain: ", basename(mcmc_file))
  fit <- readRDS(mcmc_file)
} else {
  # NOT pruned either: Figure 5a shows which of the fixed references the cohort
  # does NOT support, so the parked ones have to survive into the solution.
  fit <- SignaturePPF(
    prune_solution = FALSE,
    data,
    sigs       = CosmicSigs,      # the catalogue, held fixed
    sigs_fixed = TRUE,
    method     = "mcmc",
    prior      = prior,
    controls   = SignaturePPF_control(nsamples = NSAMPLES, burnin = BURNIN,
                                      sampler = "agess",
                                      logpost_every = LOGPOST_EVERY),
    # Initial starting point
    init = SignaturePPF_init(Theta_start  = map$Thetas,
                             Betas_start  = map$Betas,
                             Mu_start     = map$Mu,
                             Sigma2_start = map$Sigma2),
    init_mcmc_from_map = FALSE,
    prune_after_map = FALSE,
    checkpoint = SignaturePPF_checkpoint(dir = file.path(DIR_REFIT, "tmpMCMC"),
                                         every = 100, resume = TRUE),
    seed    = SEED,
    verbose = TRUE)

  saveRDS(fit, mcmc_file, compress = "gzip")
}

print(fit)

################################################################################
# 4. Posterior summaries
################################################################################
results <- posterior_summaries(fit)
saveRDS(results, file.path(DIR_REFIT, "resultsMCMC_refit.rds.gzip"),
        compress = "gzip")

mu_tbl <- data.frame(
  signature = SIGS_TO_USE,
  mu = as.numeric(fit$Mu[SIGS_TO_USE]),
  mu_low = as.numeric(results$Mu$lowCI[SIGS_TO_USE]),
  mu_high = as.numeric(results$Mu$highCI[SIGS_TO_USE]))
mu_tbl$supported <- mu_tbl$mu > 10 * fit$prior$epsilon
mu_tbl <- mu_tbl[order(-mu_tbl$mu), ]
write.csv(mu_tbl, file.path(DIR_REFIT, "signature_summary.csv"), row.names = FALSE)
print(mu_tbl)

ggsave(file.path(FIG_DIR, "Refit_logposterior_trace.pdf"),
       plot_logposterior_trace(fit), width = 7, height = 3)

message("\ndone: outputs in ", DIR_REFIT)
