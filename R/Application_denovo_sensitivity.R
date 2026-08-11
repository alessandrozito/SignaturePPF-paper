################################################################################
# Sensitivity of the de novo solution to Kmax and to the shrinkage priors
#
# Reviewer question: is the de novo solution stable when Kmax is raised beyond
# 12, and when the compressive and coefficient-shrinkage priors are varied? Do
# the same COSMIC-matched signatures, activities and covariate effects appear?
#
#   scenario     Kmax   c0    d0     what it varies
#   -----------------------------------------------------------------------
#   reference     12    100   1      the setting Application_denovo.R uses
#   K15           15    100   1      Kmax raised
#   K20           20    100   1      Kmax raised further
#   K12_c0_10     12     10   0.1    weaker beta shrinkage
#   K12_c0_1      12      2   0.02   much weaker beta shrinkage
#
# MAP ONLY. The question is whether the mode moves, and the answer does not need
# a posterior - which at 2 kb would cost a day per scenario. Each scenario gets
# the SAME treatment as the reference, three random starts with the best log
# posterior kept, so a difference between scenarios reflects the prior rather
# than the luck of one starting point.
#
# Usage:  Rscript R/Application_denovo_sensitivity.R [cores]
#
#   cores  how many scenarios to fit at once (default 3). Each worker holds the
#          2 kb cohort in memory - about 8 GB - so raise this only if the
#          machine has the room.
#
# Every start is cached to its own file, so the run is restartable and a
# scenario can be added without redoing the others.
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(parallel)
  library(ggplot2)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()
check_inputs(PATH_ICGC2KB)

args <- commandArgs(trailingOnly = TRUE)
N_CORES <- if (length(args)) as.integer(args[1]) else 3L

N_STARTS <- 3
MAXITER <- 5000
TOL <- 1e-7

################################################################################
# 1. Data
################################################################################
data <- readRDS(PATH_ICGC2KB)
v <- SignaturePPF_validate(data)
message(sprintf("%s mutations | %d samples | %d covariates | %s bins",
                format(v$N, big.mark = ","), v$J, v$p,
                format(v$nbins, big.mark = ",")))

################################################################################
# 2. The scenarios
################################################################################
# c0 and d0 are the shape and rate of the inverse-Gamma on sigma^2_k, the prior
# variance of the coefficients. Lowering c0 while lowering d0 in step keeps the
# prior mean of sigma^2 near 0.01 but fattens its tail, so the coefficients are
# progressively freer to be large.
scenarios <- list(
  reference   = list(K = K_DENOVO, c0 = 100, d0 = 1),
  K15         = list(K = 15,       c0 = 100, d0 = 1),
  K20         = list(K = 20,       c0 = 100, d0 = 1),
  K12_c0_10   = list(K = K_DENOVO, c0 = 10,  d0 = 0.1),
  K12_c0_1    = list(K = K_DENOVO, c0 = 2,   d0 = 0.02))

message("\n", length(scenarios), " scenarios x ", N_STARTS, " starts on ",
        N_CORES, " cores")

################################################################################
# 3. Fit
#
#    The reference scenario shares its fits with Application_denovo.R when that
#    has already been run: same K, same prior, same seeds, so re-optimising it
#    here would spend hours reproducing a file that is already on disk.
################################################################################
scenario_dir <- function(nm) {
  if (nm == "reference") file.path(DIR_DENOVO, "MapSolutions")
  else file.path(DIR_SENSITIVITY, nm)
}

fit_scenario <- function(nm) {
  sc <- scenarios[[nm]]
  message("\n######## scenario ", nm, " | K = ", sc$K,
          ", c0 = ", sc$c0, ", d0 = ", sc$d0, " ########")
  fit_map_restarts(
    data, out_dir = scenario_dir(nm),
    n_starts = N_STARTS, seed = SEED,
    sigs = NULL, K = sc$K,
    prior = SignaturePPF_prior(c0 = sc$c0, d0 = sc$d0),
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL),
    verbose = TRUE)
}

fits <- mclapply(names(scenarios), function(nm) try(fit_scenario(nm), silent = TRUE),
                 mc.cores = N_CORES, mc.preschedule = FALSE)
names(fits) <- names(scenarios)

failed <- names(fits)[vapply(fits, inherits, logical(1), "try-error")]
for (nm in failed) {
  message("!! scenario ", nm, " failed: ",
          conditionMessage(attr(fits[[nm]], "condition")))
}
fits <- fits[setdiff(names(fits), failed)]
stopifnot(length(fits) > 0, "reference" %in% names(fits))

for (nm in names(fits)) {
  saveRDS(fits[[nm]], file.path(DIR_SENSITIVITY, paste0("MAP_", nm, ".rds.gzip")),
          compress = "gzip")
}

################################################################################
# 4. Compare each scenario with the reference
#
#    Signatures are matched one to one against the reference by cosine, so a
#    scenario that splits or reorders them is still comparable. Reported per
#    matched pair: the cosine to the reference signature, whether both match the
#    same COSMIC signature, the correlation of the per-patient activities, and
#    the agreement of the covariate coefficients.
################################################################################
ref <- fits[["reference"]]
ref_cosmic <- match_to_cosmic(ref$Signatures)
active <- function(f) f$Mu > 10 * f$prior$epsilon

# Activities as a share of each patient's total, so scenarios with different K
# are on one scale.
share <- function(f) sweep(f$Thetas, 2, colSums(f$Thetas), "/")
ref_share <- share(ref)

compare_one <- function(nm) {
  f <- fits[[nm]]
  keep_f <- which(active(f))
  keep_r <- which(active(ref))
  m <- hungarian_match_signatures(ref$Signatures[, keep_r, drop = FALSE],
                                  f$Signatures[, keep_f, drop = FALSE])
  cos_f <- match_to_cosmic(f$Signatures)

  do.call(rbind, lapply(seq_len(nrow(m)), function(i) {
    sr <- m$signature[i]      # reference signature
    sf <- m$ppf_sig[i]        # this scenario's signature
    data.frame(
      scenario = nm,
      sig_reference = sr,
      sig_scenario = sf,
      cosine_to_reference = m$cosine[i],
      cosmic_reference = ref_cosmic$best_match[ref_cosmic$signature == sr],
      cosmic_scenario = cos_f$best_match[cos_f$signature == sf],
      activity_cor = cor(ref_share[sr, ], share(f)[sf, ]),
      beta_cor = cor(ref$Betas[, sr], f$Betas[, sf]),
      beta_rmse = sqrt(mean((ref$Betas[, sr] - f$Betas[, sf])^2)),
      beta_sign_agreement = mean(sign(ref$Betas[, sr]) == sign(f$Betas[, sf])),
      stringsAsFactors = FALSE)
  }))
}

comparison <- do.call(rbind, lapply(setdiff(names(fits), "reference"), compare_one))
comparison$same_cosmic <- comparison$cosmic_reference == comparison$cosmic_scenario
write.csv(comparison, file.path(DIR_SENSITIVITY, "sensitivity_comparison.csv"),
          row.names = FALSE)

overview <- data.frame(
  scenario = names(fits),
  K = vapply(fits, function(f) f$K, integer(1)),
  c0 = vapply(names(fits), function(n) scenarios[[n]]$c0, numeric(1)),
  d0 = vapply(names(fits), function(n) scenarios[[n]]$d0, numeric(1)),
  n_active = vapply(fits, function(f) sum(active(f)), integer(1)),
  logposterior = vapply(fits, map_logposterior, numeric(1)),
  row.names = NULL)
overview$matched_to_reference <- overview$scenario |>
  vapply(function(n) if (n == "reference") NA_integer_ else
    sum(comparison$scenario == n), integer(1))
overview$same_cosmic <- overview$scenario |>
  vapply(function(n) if (n == "reference") NA_real_ else
    mean(comparison$same_cosmic[comparison$scenario == n]), numeric(1))

write.csv(overview, file.path(DIR_SENSITIVITY, "sensitivity_overview.csv"),
          row.names = FALSE)
print(overview)
print(comparison)

saveRDS(list(overview = overview, comparison = comparison,
             scenarios = scenarios, fits = fits),
        file.path(DIR_SENSITIVITY, "sensitivity_results.rds.gzip"),
        compress = "gzip")

message("\ndone: outputs in ", DIR_SENSITIVITY)
