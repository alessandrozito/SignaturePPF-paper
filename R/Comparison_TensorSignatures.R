################################################################################
# Produces: Figure S11
#
# Comparison against TensorSignatures (Vohringer et al., Nat Commun 2021),
#
# Usage:  Rscript R/Comparison_TensorSignatures.R [rank]
#
#   rank   optional. Which TensorSignatures rank to compare against. Default is
#          the lowest-AIC rank in the sweep.
#
# Every expensive step is cached, so rerunning only redoes what is missing:
# the chromatin dataset, the SignaturePPF fit and each TensorSignatures rank are
# all skipped if their output is already on disk.
#
# PREREQUISITE - the Python environment
# -------------------------------------
# TensorSignatures 0.5.0 pins tensorflow <= 1.15, whose wheels stop at Python
# 3.7, so it cannot share an interpreter with anything modern and needs its own
# conda environment. Build it once with
#
#     ./setup_tensorsig_env.sh
#
# This script does NOT build it automatically: it downloads and installs a
# miniconda distribution under $HOME. It stops with that command if the environment is
# missing.
#
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()
check_inputs()

REFERENCE_STATE <- "Quies"
TAG <- "icgc_chromatin"
RANKS <- 4:12                 # list of ranks tested by TensorSignatures
K_PPF <- 12                   # upper bound to the number of signatures in PPF
MU_MIN <- 0.05                # PPF signatures below this are compressed, not fitted

TS_BASE <- file.path(DIR_TENSORSIG, TAG)
PATH_DATASET <- file.path(DIR_TENSORSIG, "dataset_chromatin.rds.gzip")
PATH_PPF_FIT <- file.path(DIR_TENSORSIG, "fit_ppf_chromatin.rds.gzip")

args <- commandArgs(trailingOnly = TRUE)
rank_requested <- if (length(args)) as.integer(args[1]) else NA_integer_


################################################################################
# 1. The chromatin-state dataset
################################################################################
message("\n== 1. chromatin-state dataset ==")
if (file.exists(PATH_DATASET)) {
  message("cached: ", basename(PATH_DATASET))
  dataChrom <- readRDS(PATH_DATASET)
} else {
  gr_tumor <- readRDS(PATH_ICGC_SNV)

  # Blacklisted mutations are dropped outright rather than left to bin weights:
  # a mutation inside a blacklisted region has no usable exposure behind it.
  blacklist <- rtracklayer::import(PATH_BLACKLIST)
  gr_tumor <- gr_tumor[-S4Vectors::queryHits(
    GenomicRanges::findOverlaps(gr_tumor, blacklist))]
  gr_tumor <- gr_tumor[GenomicRanges::seqnames(gr_tumor) != "chrY"]

  df_copy <- readr::read_tsv(PATH_ICGC_CN, show_col_types = FALSE)
  gr_copy <- GenomicRanges::GRanges(
    seqnames = paste0("chr", df_copy$chr),
    ranges = IRanges::IRanges(start = df_copy$start, end = df_copy$end),
    strand = "*", sample = df_copy$sampleID, score = df_copy$value)

  dataChrom <- build_chromatin_dataset(gr_tumor, gr_copy,
                                       reference = REFERENCE_STATE)
  saveRDS(dataChrom, PATH_DATASET, compress = "gzip")
}

invisible(SignaturePPF_validate(dataChrom))
message("bins: ", nrow(dataChrom$SignalTrack),
        " | states: ", ncol(dataChrom$SignalTrack) + 1,
        " | samples: ", ncol(dataChrom$CopyTrack),
        " | mutations: ", length(dataChrom$gr_Mutations))


# Count numbers of mutation by chromatin state
table(dataChrom$state_of_bin[dataChrom$bin_of_mut])
prop.table(table(dataChrom$state_of_bin[dataChrom$bin_of_mut]))

################################################################################
# 2. Fit SignaturePPF
################################################################################
message("\n== 2. SignaturePPF fit ==")
if (file.exists(PATH_PPF_FIT)) {
  message("cached: ", basename(PATH_PPF_FIT))
  fit <- readRDS(PATH_PPF_FIT)
} else {
  # NOT pruned here. This script prunes explicitly at MU_MIN further down, and
  # relabel_by_mu() runs on whatever comes back - so letting the package prune at
  # its own looser threshold would make the LABELS depend on whether the fit came
  # off disk or out of the optimizer.
  fit <- SignaturePPF(dataChrom,
                      K = K_PPF,
                      method = "map",
                      prune_solution = FALSE,
                      controls = SignaturePPF_control(maxiter = 4000, tol = 1e-6),
                      seed = SEED,
                      verbose = TRUE)
  saveRDS(fit, PATH_PPF_FIT, compress = "gzip")
}

# Relabel the signatures based on relevance weights
fit <- relabel_by_mu(fit)
print(attr(fit, "relabel"))

print(fit)
plot(fit)

ref_match <- match_to_cosmic(fit$Signatures)
ref_match$mu <- as.numeric(fit$Mu[ref_match$signature])
ref_match <- ref_match[order(-ref_match$mu), ]
write.csv(ref_match, file.path(DIR_TENSORSIG, "ppf_signature_cosmic_match.csv"),
          row.names = FALSE)
print(ref_match)

plot_Signatures(fit$Signatures[, -c(11:12)]) +  plot_Betas(x = fit$Betas[, -c(11:12)])


################################################################################
# 3. Export the same data as a TensorSignatures tensor
################################################################################
message("\n== 3. export tensor ==")
export_ts_chromatin(dataChrom, out_dir = DIR_TENSORSIG, tag = TAG)

################################################################################
# 4. Fit TensorSignatures
#
#    TensorSignatures has no automatic selection for the priors, so its number of
#    signatures has to be chosen by an explicit sweep and an information criterion.
################################################################################
message("\n== 4. TensorSignatures rank sweep ==")
if (!file.exists(TENSORSIG_PYTHON)) {
  stop("the TensorSignatures environment is missing:\n  ", TENSORSIG_PYTHON,
       "\n\nBuild it once with\n  ./setup_tensorsig_env.sh\n",
       "\nIt installs miniconda under $HOME (removable with rm -rf ~/miniconda3).",
       call. = FALSE)
}

sweep_script <- file.path(PAPER_ROOT, "run_ts_sweep.sh")
status <- system2(sweep_script, args = as.character(RANKS),
                  env = c(paste0("SIGNATUREPPF_PAPER=", shQuote(PAPER_ROOT)),
                          paste0("TENSORSIG_PYTHON=", shQuote(TENSORSIG_PYTHON)),
                          paste0("TS_TAG=", shQuote(TAG))))
if (status != 0) {
  stop("the TensorSignatures sweep exited with status ", status,
       "\n  Rerun it directly to see the full log:  ./run_ts_sweep.sh",
       call. = FALSE)
}

################################################################################
# 5. Rank selection
################################################################################
message("\n== 5. rank selection ==")
sweep <- ts_sweep_summary(TS_BASE)
if (is.null(sweep)) stop("no TensorSignatures fits under ", TS_BASE)
write.csv(sweep, file.path(DIR_TENSORSIG, "ts_rank_sweep.csv"), row.names = FALSE)
print(sweep)

# Select the model with lowest AIC
rank <- if (!is.na(rank_requested)) rank_requested else sweep$rank[which.min(sweep$AIC)]
message("using rank ", rank,
        if (is.na(rank_requested)) " (lowest AIC)" else " (given on the command line)")

TS_DIR <- file.path(TS_BASE, sprintf("rank%02d", rank))
if (!dir.exists(TS_DIR)) stop("no fit at ", TS_DIR)

sweep_long <- rbind(
  data.frame(rank = sweep$rank, criterion = "AIC", value = sweep$AIC),
  data.frame(rank = sweep$rank, criterion = "BIC", value = sweep$BIC))
chosen <- data.frame(rank = rank, criterion = "AIC",
                     value = sweep$AIC[match(rank, sweep$rank)])

p_sweep <- ggplot(sweep_long, aes(rank, value)) +
  geom_line(colour = "grey50") + geom_point() +
  geom_point(data = chosen, colour = "#CD2626", size = 3) +
  facet_wrap(~ criterion, scales = "free_y") +
  labs(x = "Rank (number of signatures)", y = NULL,
       title = "TensorSignatures rank selection",
       subtitle = paste0("selected rank ", rank, " (lowest AIC)")) +
  theme_bw()
p_sweep

################################################################################
# 6. Compare: signature spectra
################################################################################
message("\n== 6. spectra ==")
ts <- read_ts_fit(TS_DIR)

matched   <- match_and_relabel_ts(ts, fit, mu_min = MU_MIN)
ts        <- matched$ts
match_tbl <- matched$match
message("matched ", nrow(match_tbl), " pair(s) | unmatched TS: ",
        paste(attr(match_tbl, "unmatched_ts"), collapse = ", "),
        " | unmatched PPF: ",
        paste(attr(match_tbl, "unmatched_ppf"), collapse = ", "))
write.csv(match_tbl, file.path(DIR_TENSORSIG, "signature_matching.csv"),
          row.names = FALSE)
print(match_tbl)

# Each method against COSMIC, the neutral reference for "did it find a known
# signature".
cosmic_cmp <- rbind(
  cbind(method = "SignaturePPF", match_to_cosmic(fit$Signatures)),
  cbind(method = "TensorSignatures", match_to_cosmic(ts$signatures)))

################################################################################
# 7. Compare: chromatin-state effects
################################################################################

message("\n== 7. chromatin-state effects ==")

# Drop the signatures the compressive prior switched off. Same threshold the
# matching above used, so the panels and the matching describe the same set.
fit_filter <- prune_signatures(fit, threshold = MU_MIN)

# Plot Signatures side by side, and beta coefficients for PPF
p_sig_ts <- (plot_Signatures(matched$ts$signatures) +
               theme(axis.text.x = element_blank()))
p_sig_PPF <- plot_Signatures(fit_filter$Signatures) +
  theme(axis.text.x = element_blank())
p_betas_PPF <- plot_Betas(x = fit_filter$Betas) +
  theme(plot.margin = margin(r = 20, unit = "pt"))

# Filter fit object for plots
p_mu <- plot_Mu(fit_filter, dataChrom) + theme(legend.position = "right")

ggsave(file.path(FIG_DIR, "TensorSignatures_top_panel_mu.pdf"),
       p_mu,
       width = 3.23, height = 4.33)


# Display and save the three plots displayed
p_all_plots <- p_sig_ts + p_sig_PPF + p_betas_PPF + plot_layout(widths = c(1,1,2))
p_all_plots
ggsave(file.path(FIG_DIR, "TensorSignatures_top_panel.pdf"),
       p_all_plots,
       width = 10.82, height = 5.84)

# Match to COSMIC
match_to_cosmic(fit_filter$Signatures)
cosine(fit_filter$Signatures, SignaturePPF::COSMIC_v3.4_SBS96_GRCh37[, c("SBS36", "SBS18")])

# Plot betas vs TS amplitudes
cmp <- compare_chromatin_effects(fit, ts, match_tbl, reference = REFERENCE_STATE)
write.csv(cmp, file.path(DIR_TENSORSIG, "chromatin_effect_comparison.csv"),
          row.names = FALSE)

p_by_signature <- plot_chromatin_effects(cmp, by = "signature", ncol = 5)
ggsave(file.path(FIG_DIR, "TensorSignatures_chromatin_effects.pdf"),
       p_by_signature, width = 10.82, height = 3.65)

effect_agreement <- data.frame(
  n_pairs = length(unique(cmp$pair_label)),
  n_points = nrow(cmp),
  pearson = cor(cmp$beta, cmp$ts_logratio),
  spearman = cor(cmp$beta, cmp$ts_logratio, method = "spearman"),
  sign_agreement = mean(sign(cmp$beta) == sign(cmp$ts_logratio)))
print(effect_agreement)

message("\ndone: tables in ", DIR_TENSORSIG, "\n      figures in ", FIG_DIR)
