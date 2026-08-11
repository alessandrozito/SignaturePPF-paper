################################################################################
# De novo application: ICGC Breast-AdenoCa at 2 kb
#
# K = 12 signatures estimated from the data. Three random starting points for
# the MAP, keep the one with the highest log posterior, then sample from there.
#
#   MAP     3 starts, K = 12, cached one file each
#   MCMC    10000 sweeps, 5000 burn-in, started at the best MAP
#
# Signatures the compressive prior parks near epsilon are NOT pruned before
# sampling: whether a marginal signature survives is one of the things the
# posterior is being asked, and fixing the dimension at the mode's answer would
# decide it in advance.
#
# The chain is checkpointed every 100 sweeps, so an interrupted run resumes
# bit-exactly - the RNG state is saved with each block.
#
# Usage:  Rscript R/Application_denovo.R
#
# BEWARE: at 2 kb this is 1.39 million bins and 707k mutations. One sweep costs
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

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()
check_inputs(PATH_ICGC2KB)

## ------------------------------------------------------------------ settings
N_STARTS <- 3
MAXITER <- 5000
TOL <- 1e-7
NSAMPLES <- 10000
BURNIN <- 5000
LOGPOST_EVERY <- 10      # the log posterior is a diagnostic; evaluating it every
                         # sweep costs a full pass over the data for nothing

################################################################################
# 1. Data
################################################################################
data <- readRDS(PATH_ICGC2KB)
v <- SignaturePPF_validate(data)
message(sprintf("%s mutations | %d samples | %d covariates | %s bins",
                format(v$N, big.mark = ","), v$J, v$p,
                format(v$nbins, big.mark = ",")))

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
# 3. MCMC from the best MAP
################################################################################
message("\n== MCMC ==")
fit <- run_mcmc_from_map(
  data, map, out_dir = DIR_DENOVO,
  sigs = NULL, sigs_fixed = FALSE,
  K = K_DENOVO,
  prior = prior,
  controls = SignaturePPF_control(nsamples = NSAMPLES, burnin = BURNIN,
                                  logpost_every = LOGPOST_EVERY),
  every = 100, seed = SEED)

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
