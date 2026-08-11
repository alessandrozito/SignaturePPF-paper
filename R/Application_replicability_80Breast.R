################################################################################
# Replication: 80 breast cancers (Davies et al. 2017) vs ICGC Breast-AdenoCa
#
# Both cohorts are binned at 10 kb (rather than the 2 kb used for the main
# application) and carry the same 11 covariates. We only run the MAP estimate,
# keeping signatures fixed.
#
# Usage:  Rscript R/Application_replicability_80Breast.R
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()
check_inputs()

FIT_80   <- file.path(DIR_REPLICATION, "fit_Breast80_10kb_map.rds.gzip")
FIT_ICGC <- file.path(DIR_REPLICATION, "fit_BreastICGC_10kb_map.rds.gzip")

################################################################################
# 1. Data
################################################################################
data80 <- load_cohort(PATH_BREAST80)
dataICGC <- load_cohort(PATH_ICGC10KB)

# Fail here rather than three minutes into an optimisation.
invisible(SignaturePPF_validate(data80))
invisible(SignaturePPF_validate(dataICGC))

stopifnot(identical(colnames(data80$SignalTrack), colnames(dataICGC$SignalTrack)),
          nrow(data80$SignalTrack) == nrow(dataICGC$SignalTrack))

CosmicSigs <- COSMIC_v3.4_SBS96_GRCh37[, SIGS_TO_USE]

################################################################################
# 2. Fit both cohorts
#
#    We keep `sigs_fixed = TRUE`
################################################################################
controls <- SignaturePPF_control(maxiter = 200, tol = 1e-6)

fit_ppf <- function(data, out_file) {
  if (file.exists(out_file)) {
    message("using existing fit: ", basename(out_file))
    return(readRDS(out_file))
  }
  fit <- SignaturePPF(data,
                      sigs = CosmicSigs,
                      sigs_fixed = TRUE,
                      method = "map",
                      controls = controls,
                      seed = SEED,
                      verbose = TRUE)
  saveRDS(fit, out_file, compress = "gzip")
  fit
}

# Run both models and save the output
fitPPF_80 <- fit_ppf(data80, FIT_80)
fitPPF_ICGC <- fit_ppf(dataICGC, FIT_ICGC)

print(fitPPF_80)
print(fitPPF_ICGC)

################################################################################
# 3. Detect differences in regression coefficients in each cohort
################################################################################
mu_table <- data.frame(
  signature = colnames(fitPPF_80$Signatures),
  mu_Breast80 = as.numeric(fitPPF_80$Mu[colnames(fitPPF_80$Signatures)]),
  mu_ICGC = as.numeric(fitPPF_ICGC$Mu[colnames(fitPPF_80$Signatures)]),
  row.names = NULL)
mu_table$assigned_Breast80 <- df_assign(fitPPF_80, data80)$m[
  match(mu_table$signature, df_assign(fitPPF_80, data80)$best_sig)]
mu_table$assigned_ICGC <- df_assign(fitPPF_ICGC, dataICGC)$m[
  match(mu_table$signature, df_assign(fitPPF_ICGC, dataICGC)$best_sig)]

write.csv(mu_table, file.path(DIR_REPLICATION, "relevance_weights.csv"),
          row.names = FALSE)
print(mu_table)

################################################################################
# 4. Figures
################################################################################
p_burden <- plot_burden_along_genome(list(Breast80 = data80, BreastICGC = dataICGC))
ggsave(file.path(FIG_DIR, "Replication_burden_along_genome.pdf"), p_burden,
       width = 11, height = 3.2)

# --- per-cohort coefficient heatmaps, each beside its relevance-weight column.
#     cap = 1 bounds the colour scale only; the printed number is the estimate.
panel <- function(fit, data) {
  plot_Betas(fit, cap = 1) +
    plot_vector_facets_x(df_assign(fit, data), levs = colnames(fit$Betas)) +
    plot_layout(widths = c(4, 1))
}
p_80 <- panel(fitPPF_80, data80)
p_icgc <- panel(fitPPF_ICGC, dataICGC)

ggsave(file.path(FIG_DIR, "Replication_betas_Breast80.pdf"), p_80, width = 9, height = 6)
ggsave(file.path(FIG_DIR, "Replication_betas_BreastICGC.pdf"), p_icgc, width = 9, height = 6)

# --- the replication figure itself
p_rep <- plot_beta_replication(fitPPF_80, fitPPF_ICGC,
                               label_x = "80 Breast", label_y = "ICGC")
ggsave(file.path(FIG_DIR, "Replication_betas_scatter.pdf"), p_rep,
       width = 7, height = 5.5)

# --- and the difference, on the same capped scale as the two panels
p_diff <- plot_Betas(fitPPF_80$Betas - fitPPF_ICGC$Betas, cap = 1) +
  ggtitle("Breast80 - ICGC")
ggsave(file.path(FIG_DIR, "Replication_betas_difference.pdf"), p_diff, width = 7, height = 6)

################################################################################
# 5. Agreement summary
#
#    Restricted to signatures the compressive prior kept in BOTH cohorts: for a
#    switched-off signature beta is a draw from its prior, so including it would
#    measure the prior rather than the replication.
################################################################################
shared <- intersect(
  names(which(fitPPF_80$Mu > 0.01)),
  names(which(fitPPF_ICGC$Mu > 0.01)))

x <- as.numeric(fitPPF_80$Betas[, shared])
y <- as.numeric(fitPPF_ICGC$Betas[, shared])
agreement <- data.frame(
  n_signatures_shared = length(shared),
  n_coefficients = length(x),
  pearson = cor(x, y),
  spearman = cor(x, y, method = "spearman"),
  sign_agreement = mean(sign(x) == sign(y)),
  rmse = sqrt(mean((x - y)^2)))
write.csv(agreement, file.path(DIR_REPLICATION, "beta_agreement.csv"),
          row.names = FALSE)
print(agreement)

message("done: outputs in ", DIR_REPLICATION, " and ", FIG_DIR)
