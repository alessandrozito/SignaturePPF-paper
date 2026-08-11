################################################################################
# TensorSignatures comparison, step 3 of 3: import the fits and compare.
#
# Requires R/02_tensorsignatures_prepare.R and then the Python rank sweep
# (bash/run_ts_sweep.sh) to have run.
#
# Three comparisons, in increasing order of what they actually test:
#
#   1. SPECTRA   - do the two methods find the same signatures? Matched
#                  one-to-one by cosine similarity (Hungarian, not argmax).
#   2. EFFECTS   - do they agree on the chromatin-state effect of each matched
#                  signature? Both are log enrichments relative to the reference
#                  state, so they are directly comparable with no rescaling.
#   3. RATE      - do they predict WHERE the mutations are? This is where the
#                  models genuinely differ: PPF has an intensity per bin, while
#                  TensorSignatures can only place a total per (state, sample)
#                  and has no notion of position within a state.
#
# Usage:  Rscript R/03_tensorsignatures_compare.R [rank]
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()

REFERENCE_STATE <- "Quies"
TAG <- "icgc_chromatin"
TS_BASE <- file.path(DIR_TENSORSIG, TAG)

dat <- readRDS(file.path(DIR_TENSORSIG, "dataset_chromatin.rds.gzip"))
fit <- readRDS(file.path(DIR_TENSORSIG, "fit_ppf_chromatin.rds.gzip"))

################################################################################
# 1. Rank selection
#
#    TensorSignatures has no compressive prior, so its number of signatures is
#    chosen by information criterion over an explicit sweep. PPF selects K
#    through the prior instead, which is one of the differences worth reporting.
################################################################################
sweep <- ts_sweep_summary(TS_BASE)
if (is.null(sweep)) {
  stop("no TensorSignatures fits under ", TS_BASE,
       "\n  Run bash/run_ts_sweep.sh first.")
}
write.csv(sweep, file.path(DIR_TENSORSIG, "ts_rank_sweep.csv"), row.names = FALSE)
print(sweep)

args <- commandArgs(trailingOnly = TRUE)
rank <- if (length(args)) as.integer(args[1]) else sweep$rank[which.min(sweep$BIC)]
message("using TensorSignatures rank ", rank,
        if (!length(args)) " (lowest BIC)" else " (given on the command line)")

TS_DIR <- file.path(TS_BASE, sprintf("rank%02d", rank))
if (!dir.exists(TS_DIR)) stop("no fit at ", TS_DIR)

p_sweep <- ggplot(sweep, aes(rank, BIC)) +
  geom_line(colour = "grey50") + geom_point() +
  geom_point(data = sweep[which.min(sweep$BIC), ], colour = "#CD2626", size = 3) +
  labs(x = "Rank (number of signatures)", y = "BIC",
       title = "TensorSignatures rank selection") +
  theme_bw()
ggsave(file.path(FIG_DIR, "03_ts_rank_sweep.pdf"), p_sweep, width = 5, height = 3.5)

################################################################################
# 2. Spectra
################################################################################
ts <- read_ts_fit(TS_DIR)
match_tbl <- hungarian_match_signatures(ts$signatures, fit$Signatures)
message("matched ", nrow(match_tbl), " signature pair(s); ",
        "unmatched TS: ", paste(attr(match_tbl, "unmatched_ts"), collapse = ", "),
        " | unmatched PPF: ",
        paste(attr(match_tbl, "unmatched_ppf"), collapse = ", "))

# Each method against COSMIC, which is the neutral reference for "did it find a
# known signature".
cosmic_cmp <- rbind(
  cbind(method = "SignaturePPF", match_to_cosmic(fit$Signatures)),
  cbind(method = "TensorSignatures", match_to_cosmic(ts$signatures)))
write.csv(cosmic_cmp, file.path(DIR_TENSORSIG, "signature_cosmic_comparison.csv"),
          row.names = FALSE)
write.csv(match_tbl, file.path(DIR_TENSORSIG, "signature_matching.csv"),
          row.names = FALSE)
print(cosmic_cmp)

p_cos <- ggplot(cosmic_cmp, aes(method, cosine)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92") +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.8) +
  labs(x = NULL, y = "Best cosine similarity to COSMIC v3.4") +
  theme_bw()
ggsave(file.path(FIG_DIR, "03_cosine_to_cosmic.pdf"), p_cos, width = 4.5, height = 4)

################################################################################
# 3. Chromatin-state effects
################################################################################
cmp <- compare_chromatin_effects(fit, TS_DIR, reference = REFERENCE_STATE)
write.csv(cmp, file.path(DIR_TENSORSIG, "chromatin_effect_comparison.csv"),
          row.names = FALSE)

ggsave(file.path(FIG_DIR, "03_chromatin_effects_pooled.pdf"),
       plot_chromatin_effects(cmp), width = 6, height = 5)
ggsave(file.path(FIG_DIR, "03_chromatin_effects_by_signature.pdf"),
       plot_chromatin_effects(cmp, by_signature = TRUE), width = 12, height = 9)

effect_agreement <- data.frame(
  n_pairs = length(unique(cmp$pair_label)),
  n_points = nrow(cmp),
  pearson = cor(cmp$beta, cmp$ts_logratio),
  spearman = cor(cmp$beta, cmp$ts_logratio, method = "spearman"),
  sign_agreement = mean(sign(cmp$beta) == sign(cmp$ts_logratio)))
write.csv(effect_agreement,
          file.path(DIR_TENSORSIG, "chromatin_effect_agreement.csv"),
          row.names = FALSE)
print(effect_agreement)

################################################################################
# 4. Regional mutation rate
#
#    The comparison that separates the two models: aggregated into 1 Mb windows,
#    how close is each method's predicted burden to the observed one?
################################################################################
rate <- compare_mutation_rate(dat, fit, TS_DIR, window = 1e6)
write.csv(rate, file.path(DIR_TENSORSIG, "mutation_rate_windows.csv"),
          row.names = FALSE)

scores <- score_mutation_rate(rate)
write.csv(scores, file.path(DIR_TENSORSIG, "mutation_rate_scores.csv"),
          row.names = FALSE)
print(scores)

ggsave(file.path(FIG_DIR, "03_mutation_rate_along_genome.pdf"),
       plot_mutation_rate(rate), width = 12, height = 4)

message("done: outputs in ", DIR_TENSORSIG, " and ", FIG_DIR)
