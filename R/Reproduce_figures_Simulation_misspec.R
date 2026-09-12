################################################################################
# Produces: Figures S2 and S3
#
# Figures for the misspecification simulation study
#
# Reads output/Simulation_misspec/all_calibration.csv (written by the `score`
# stage of Simulation_misspec.R) and draws the expected calibration error of the
# mutation-level attribution probabilities, in and out of sample.
#
#   Rscript R/Reproduce_figures_Simulation_misspec.R
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))

calib <- read.csv(file.path(DIR_SIM_MISSPEC, "all_calibration.csv"))

## ------------------------------------------------------------------- labels
MODEL_LABELS <- c(CompNMF_Fixed = "CompressiveNMF",
                  map_Fixed     = "PPF (MAP)",
                  mcmc_Fixed    = "PPF (MCMC)")
MODEL_COLS <- c("CompressiveNMF" = "#CD2626",
                "PPF (MAP)"      = "#FF8C00",
                "PPF (MCMC)"     = "#000D8B")

SCENARIO_LABELS <- c(S0_baseline       = "Scenario 0",
                     S1_epigenome_v025 = "Scenario 1",
                     S2_epigenome_v1   = "Scenario 2",
                     S3_epigenome_v4   = "Scenario 3",
                     S4_hotspots       = "Scenario 4",
                     S5_cn_noise       = "Scenario 5",
                     S6_opportunity    = "Scenario 6")

calib <- calib %>%
  mutate(Model = factor(MODEL_LABELS[.data$model], levels = unname(MODEL_LABELS)),
         Scenario = factor(SCENARIO_LABELS[.data$Scenario],
                           levels = unname(SCENARIO_LABELS)))

## ------------------------------------------------------- de novo results
results <- read.csv(file.path(DIR_SIM_MISSPEC, "all_results.csv"))

DENOVO_LABELS <- c(CompNMF           = "CompressiveNMF",
                   SignatureAnalyzer = "SignatureAnalyzer",
                   PPF_map           = "PPF (MAP)",
                   PPF_mcmc          = "PPF (MCMC)")
DENOVO_COLS <- c("CompressiveNMF"    = "#CD2626",
                 "SignatureAnalyzer" = "antiquewhite3",
                 "PPF (MAP)"         = "#FF8C00",
                 "PPF (MCMC)"        = "#000D8B")

# The covariate-free competitors carry no B, so their "rmse_Betas" is just the
# norm of the true effects. Blanked rather than drawn, as in Simulation_main.R.
NO_BETAS <- c("CompNMF", "SignatureAnalyzer")

results <- results %>%
  mutate(Model = factor(DENOVO_LABELS[.data$model], levels = unname(DENOVO_LABELS)),
         Scenario = factor(SCENARIO_LABELS[.data$Scenario],
                           levels = unname(SCENARIO_LABELS)))


## ------------------------------- Figure: parameter and intensity recovery
recovery <- results %>%
  mutate(K = .data$Kest,
         F1 = .data$F1,
         `Sign recovery of B` = ifelse(.data$model %in% NO_BETAS, NA,
                                       .data$sign_Betas),
         Signatures = .data$rmse_sig,
         `RMSE counts, in sample` = .data$rmse_lambda_in,
         `RMSE counts, out of sample` = .data$rmse_lambda_out,
         `Theta` = .data$rmse_theta,
         Betas = ifelse(.data$model %in% NO_BETAS, NA, .data$rmse_Betas)) %>%
  dplyr::select("Model", "K", "F1", "Sign recovery of B", "Scenario", "Simulation", "Signatures",
                "RMSE counts, in sample", "RMSE counts, out of sample",
                "Theta", "Betas") %>%
  tidyr::gather("key", "value", -"Model", -"Scenario", -"Simulation") %>%
  mutate(key = factor(.data$key,
                      levels = c("K", "F1", "Signatures", "Theta", "Betas",
                                 "Sign recovery of B",
                                 "RMSE counts, in sample",
                                 "RMSE counts, out of sample")))

p_recovery <- ggplot(recovery) +
  geom_boxplot(aes(x = .data$Scenario, y = .data$value,
                   fill = .data$Model, colour = .data$Model),
               alpha = 0.6, outlier.size = 0.4, linewidth = 0.3) +
  facet_wrap(~ key, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  scale_fill_manual(name = "Model", values = DENOVO_COLS) +
  scale_colour_manual(name = "Model", values = DENOVO_COLS) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Scenario")

ggsave(file.path(FIG_DIR, "Simulation_misspec_recovery.pdf"), p_recovery,
       width = 10.85, height = 5.10)


## ------------------------ Figure: attribution accuracy and ECE, in and out
cal_long <- calib %>%
  dplyr::select("Model", "Scenario", "Simulation",
                `Accuracy, in sample`     = "attribution_acc_in",
                `Accuracy, out of sample` = "attribution_acc_out",
                `ECE, in sample`          = "ece_in",
                `ECE, out of sample`      = "ece_out") %>%
  tidyr::gather("key", "value", -"Model", -"Scenario", -"Simulation") %>%
  mutate(key = factor(.data$key,
                      levels = c("Accuracy, in sample", "Accuracy, out of sample",
                                 "ECE, in sample", "ECE, out of sample")))

# The accuracy is a proportion and the ECE spans an order of magnitude, so the
# two pairs of panels get different y scales rather than one shared transform.
p_ece <- ggplot(cal_long) +
  geom_boxplot(aes(x = .data$Scenario, y = .data$value,
                   fill = .data$Model, colour = .data$Model),
               alpha = 0.6, outlier.size = 0.5, linewidth = 0.3) +
  facet_wrap(~ key, nrow = 1, scales = "free_y") +
  ggh4x::facetted_pos_scales(y = list(
    scale_y_continuous(),
    scale_y_continuous(),
    scale_y_log10(),
    scale_y_log10())) +
  scale_fill_manual(name = "Model", values = MODEL_COLS) +
  scale_colour_manual(name = "Model", values = MODEL_COLS) +
  theme_bw() +
  theme(aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Scenario")

ggsave(file.path(FIG_DIR, "Simulation_misspec_ECE.pdf"), p_ece,
       width = 11.4, height = 3.0)

## ---- Cumulative calibration curves: one dataset per scenario, PPF vs NMF
# One replicate per scenario, both refitted on the fixed catalogue: PPF on the
# top row, the covariate-free NMF baseline on the bottom, 2 x 7 square panels.
CAL_REPLICATE <- "Simulation_01"
CAL_SET <- "all"                # "all", "in" or "out"
CAL_THIN <- 2000                # points drawn per curve
CAL_FITS <- c("PPF" = "output_mcmc_Fixed.rds.gzip",
              "Baseline NMF" = "output_CompNMF_Fixed.rds.gzip")

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
})
source(file.path(R_DIR, "Simulation_functions.R"))
source(file.path(R_DIR, "Simulation_functions_misspec.R"))

curves <- do.call(rbind, lapply(names(SCENARIO_LABELS), function(sc) {
  dir_sc <- file.path(DIR_SIM_MISSPEC, sc, CAL_REPLICATE)
  dat <- readRDS(file.path(dir_sc, "data.rds.gzip"))
  do.call(rbind, lapply(names(CAL_FITS), function(mod) {
    fit <- open_rds_file(file.path(dir_sc, CAL_FITS[[mod]]))
    if (is.null(fit)) return(NULL)
    cbind(calibration_curve_data(dat, fit, set = CAL_SET, thin = CAL_THIN),
          Scenario = SCENARIO_LABELS[[sc]], Model = mod)
  }))
})) %>%
  mutate(Scenario = factor(.data$Scenario, levels = unname(SCENARIO_LABELS)),
         Model = factor(.data$Model, levels = names(CAL_FITS)))

curve_ends <- calibration_curve_ends(curves, by = c("type", "Scenario", "Model"))

p_curves <- ggplot(curves, aes(x = .data$conf, y = .data$pcorrect,
                               colour = .data$type)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50",
              linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(data = curve_ends, size = 1.1) +
  facet_grid(Model ~ Scenario) +
  scale_colour_manual(values = CALIBRATION_COLS, name = "Mutation type") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.spacing.x = grid::unit(0.5, "lines"),
        axis.text = element_text(size = 6.5)) +
  labs(x = "Cumulative mean confidence", y = "Cumulative P(correct)")

ggsave(file.path(FIG_DIR, "Simulation_misspec_calibration_curves.pdf"), p_curves,
       width = 11.4, height = 3.7)

p_ece / p_curves
ggsave(file.path(FIG_DIR, "Simulation_misspec_calibration_curves_ECE.pdf"), p_ece / p_curves,
       width = 11.65, height = 6.89)


## ------------------------------------------------------------ summary table
summarise_metric <- function(x) sprintf("%.4f (%.4f, %.4f)",
                                        median(x), quantile(x, 0.25), quantile(x, 0.75))

tab <- calib %>%
  group_by(.data$Scenario, .data$Model) %>%
  summarise(ECE_in = summarise_metric(.data$ece_in),
            ECE_out = summarise_metric(.data$ece_out),
            Brier_in = summarise_metric(.data$brier_in),
            Brier_out = summarise_metric(.data$brier_out),
            Accuracy_in = summarise_metric(.data$attribution_acc_in),
            Accuracy_out = summarise_metric(.data$attribution_acc_out),
            .groups = "drop")

write.csv(tab, file.path(DIR_SIM_MISSPEC, "table_calibration_summary.csv"),
          row.names = FALSE)

tab_denovo <- results %>%
  group_by(.data$Scenario, .data$Model) %>%
  summarise(Kest = median(.data$Kest),
            F1 = summarise_metric(.data$F1),
            cosine = summarise_metric(.data$cosine_R),
            rmse_sig = summarise_metric(.data$rmse_sig),
            rmse_theta = summarise_metric(.data$rmse_theta),
            rmse_Betas = ifelse(unique(.data$Model) %in%
                                  DENOVO_LABELS[NO_BETAS], NA,
                                summarise_metric(.data$rmse_Betas)),
            rmse_lambda_in = summarise_metric(.data$rmse_lambda_in),
            rmse_lambda_out = summarise_metric(.data$rmse_lambda_out),
            .groups = "drop")

write.csv(tab_denovo, file.path(DIR_SIM_MISSPEC, "table_denovo_summary.csv"),
          row.names = FALSE)

message("figures written to ", FIG_DIR)
print(as.data.frame(tab[, c("Scenario", "Model", "ECE_in", "ECE_out")]), row.names = FALSE)
