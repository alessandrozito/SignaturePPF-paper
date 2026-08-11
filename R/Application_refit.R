################################################################################
# Refit application: ICGC Breast-AdenoCa at 2 kb, signatures held fixed
#
# The fifteen COSMIC signatures with support in breast cohorts are held at their
# catalogue values and only the activities and the covariate coefficients are
# estimated.
#
#   MAP     one run - with the signatures fixed the objective is far better
#           behaved than the de novo one, so multiple starts buy nothing
#   MCMC    10000 sweeps, 5000 burn-in, started at the MAP
#
# The compressive prior still applies, so a catalogue signature the cohort does
# not support is parked near epsilon. Those are NOT pruned before sampling:
# which of the fifteen the data actually supports is a result, not a setting.
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

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
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
v <- SignaturePPF_validate(data)
message(sprintf("%s mutations | %d samples | %d covariates | %s bins",
                format(v$N, big.mark = ","), v$J, v$p,
                format(v$nbins, big.mark = ",")))

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
  map <- SignaturePPF(data,
                      sigs = CosmicSigs, sigs_fixed = TRUE,
                      method = "map",
                      prior = prior,
                      controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
                      seed = SEED, verbose = TRUE)
  saveRDS(map, map_file, compress = "gzip")
}
print(map)

################################################################################
# 3. MCMC from the MAP
################################################################################
message("\n== MCMC ==")
fit <- run_mcmc_from_map(
  data, map, out_dir = DIR_REFIT,
  sigs = CosmicSigs, sigs_fixed = TRUE,
  prior = prior,
  controls = SignaturePPF_control(nsamples = NSAMPLES, burnin = BURNIN,
                                  logpost_every = LOGPOST_EVERY),
  every = 100, seed = SEED)

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
