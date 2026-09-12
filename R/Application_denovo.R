################################################################################
# Produces: no figure. Fits the de novo model, which
#           Figures 2, 3, S5, S8 and S12 are drawn from
#
# De novo application: ICGC Breast-AdenoCa at 2 kb
#
# K = 12 signatures estimated from the data. Three random starting points for
# the MAP, keep the one with the highest log posterior, then sample from there.
#
#   MAP     3 starts, K = 12, cached one file each
#   MCMC    10000 iterations, 5000 burn-in, started at the best MAP
#
# The chain is checkpointed every 100 iterations, so an interrupted run resumes
# bit-exactly - the RNG state is saved with each block.
#
# Usage:  Rscript R/Application_denovo.R
#
# BEWARE: at 2 kb this is 1.39 million bins and 707k mutations. One iteration costs
# a few seconds, so the chain is a run of many hours. Launch it detached:
#
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nohup setsid \
#     taskset -c 0 Rscript R/Application_denovo.R \
#     > output/Application_denovo/run.log 2>&1 < /dev/null &
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
N_STARTS <- 3
MAXITER <- 5000
TOL <- 1e-7
NSAMPLES <- 10000
BURNIN <- 5000
LOGPOST_EVERY <- 10      # the log posterior is a diagnostic; evaluating it every
                         # iteration costs a full pass over the data for nothing

################################################################################
# 1. Data
################################################################################
data <- readRDS(PATH_ICGC2KB)
message(sprintf("%s mutations | %d samples | %d covariates | %s bins",
                format(length(data$gr_Mutations), big.mark = ","),
                ncol(data$CopyTrack), ncol(data$SignalTrack),
                format(nrow(data$SignalTrack), big.mark = ",")))

prior <- SignaturePPF_prior()          # a = alpha = 1.01, epsilon = 0.001, c0 = 100, d0 = 1

################################################################################
# 2. MAP from three random starts
################################################################################
message("\n== MAP ==")
map <- fit_map_restarts(
  data, out_dir = file.path(DIR_DENOVO, "MapSolutions"),
  n_starts = N_STARTS, seed = SEED,
  sigs = NULL, K = K_DENOVO,
  prior = prior,
  controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL))

saveRDS(map, file.path(DIR_DENOVO, "MAPSolution.rds.gzip"), compress = "gzip")
print(map)

################################################################################
# 3. MCMC, started from the MAP with the highest logposterior
################################################################################
message("\n== MCMC ==")
message(sprintf("starting the chain from MAP start %d of %d, log posterior %s",
                map$start, N_STARTS,
                format(map_logposterior(map), big.mark = ",", nsmall = 2)))

mcmc_file <- file.path(DIR_DENOVO, "MCMCSolution.rds.gzip")

if (file.exists(mcmc_file)) {
  message("using the existing chain: ", basename(mcmc_file))
  fit <- readRDS(mcmc_file)
} else {
  fit <- SignaturePPF(
    data,
    sigs       = NULL,            # de novo: the signatures are estimated too
    sigs_fixed = FALSE,
    K          = K_DENOVO,
    method     = "mcmc",
    prior      = prior,
    controls   = SignaturePPF_control(nsamples = NSAMPLES, burnin = BURNIN,
                                      sampler = "agess",
                                      logpost_every = LOGPOST_EVERY),

    init = SignaturePPF_init(R_start      = map$Signatures,
                             Theta_start  = map$Thetas,
                             Betas_start  = map$Betas,
                             Mu_start     = map$Mu,
                             Sigma2_start = map$Sigma2),
    init_mcmc_from_map = FALSE,
    prune_after_map = FALSE,
    checkpoint = SignaturePPF_checkpoint(dir = file.path(DIR_DENOVO, "tmpMCMC"),
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
saveRDS(results, file.path(DIR_DENOVO, "resultsMCMC_denovo.rds.gzip"),
        compress = "gzip")

mu_tbl <- data.frame(
  signature = colnames(fit$Signatures),
  mu = as.numeric(fit$Mu),
  mu_low = as.numeric(results$Mu$lowCI),
  mu_high = as.numeric(results$Mu$highCI),
  best_cosmic = match_to_cosmic(fit$Signatures)$best_match,
  cosine = match_to_cosmic(fit$Signatures)$cosine)
mu_tbl <- mu_tbl[order(-mu_tbl$mu), ]
write.csv(mu_tbl, file.path(DIR_DENOVO, "signature_summary.csv"), row.names = FALSE)
print(mu_tbl)

ggsave(file.path(FIG_DIR, "Denovo_logposterior_trace.pdf"),
       plot_logposterior_trace(fit), width = 7, height = 3)

message("\ndone: outputs in ", DIR_DENOVO)
