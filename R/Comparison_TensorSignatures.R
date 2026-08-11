################################################################################
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
# -----------------------------------------------
# Both extend 96-channel NMF with genomic covariates, but parameterise the
# genomic dependence differently - continuous log-linear against discrete
# per-state amplitudes - and only TensorSignatures models strand asymmetry, only
# PPF models copy number. We make the comparison using the ChromHMM 15-state
# annotation of breast epithelium. The bins ARE the ChromHMM segments, so the
# state assignment is exact for both methods, and PPF is given the states as
# one-hot covariates with `Quies` dropped as reference - which makes each beta
# the log enrichment relative to Quies, the same quantity TensorSignatures
# reports as a state amplitude.
#
# Two comparisons:
#
#   SPECTRA  do the two methods find the same signatures? Matched one-to-one by
#            cosine similarity, using the hungarian algorithm
#   EFFECTS  do they agree on the chromatin-state effect of each matched
#            signature? Both are log enrichments against the same reference
#            state, so they compare directly with no rescaling.
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
check_inputs()

REFERENCE_STATE <- "Quies"
TAG <- "icgc_chromatin"
RANKS <- 4:12                 # list of ranks tested by TensorSignatures
K_PPF <- 12                   # upper bound to the number of signatures in PPF

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

################################################################################
# 2. Fit SignaturePPF
################################################################################
message("\n== 2. SignaturePPF fit ==")
if (file.exists(PATH_PPF_FIT)) {
  message("cached: ", basename(PATH_PPF_FIT))
  fit <- readRDS(PATH_PPF_FIT)
} else {
  fit <- SignaturePPF(dataChrom,
                      K = K_PPF,
                      method = "map",
                      controls = SignaturePPF_control(maxiter = 4000, tol = 1e-6),
                      seed = SEED,
                      verbose = TRUE)
  saveRDS(fit, PATH_PPF_FIT, compress = "gzip")
}
print(fit)
plot(fit, what = "betas") + plot(fit, what = "mu")
plot(fit)

ref_match <- match_to_cosmic(fit$Signatures)
ref_match$mu <- as.numeric(fit$Mu[ref_match$signature])
ref_match <- ref_match[order(-ref_match$mu), ]
write.csv(ref_match, file.path(DIR_TENSORSIG, "ppf_signature_cosmic_match.csv"),
          row.names = FALSE)
print(ref_match)

ggsave(file.path(FIG_DIR, "TensorSignatures_PPF_chromatin_betas.pdf"),
       plot_chromatin_betas(fit, reference = REFERENCE_STATE),
       width = 11, height = 8)

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

# AIC, not BIC. Both are reported by the Python side and both are plotted, but
# the default is AIC: the parameter count here is large (4*95 spectrum
# parameters per signature, plus one exposure per signature-sample pair) while
# `observations` counts every cell of the count tensor, most of which are zero.
# BIC's log(n) penalty is therefore severe enough to keep selecting a rank below
# the point where the fit stops improving, which understates the signature set
# TensorSignatures would actually be run with. Pass a rank explicitly to override.
rank <- if (!is.na(rank_requested)) rank_requested else sweep$rank[which.min(sweep$AIC)]
message("using rank ", rank,
        if (is.na(rank_requested)) " (lowest AIC)" else " (given on the command line)")

TS_DIR <- file.path(TS_BASE, sprintf("rank%02d", rank))
if (!dir.exists(TS_DIR)) stop("no fit at ", TS_DIR)

# Both criteria, so the choice is visible rather than asserted; the selected
# rank is marked on the panel it was chosen from.
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
ggsave(file.path(FIG_DIR, "TensorSignatures_rank_sweep.pdf"), p_sweep,
       width = 8, height = 3.5)

################################################################################
# 6. Compare: signature spectra
################################################################################
message("\n== 6. spectra ==")
ts <- read_ts_fit(TS_DIR)
match_tbl <- hungarian_match_signatures(ts$signatures, fit$Signatures)
message("matched ", nrow(match_tbl), " pair(s) | unmatched TS: ",
        paste(attr(match_tbl, "unmatched_ts"), collapse = ", "),
        " | unmatched PPF: ",
        paste(attr(match_tbl, "unmatched_ppf"), collapse = ", "))
write.csv(match_tbl, file.path(DIR_TENSORSIG, "signature_matching.csv"),
          row.names = FALSE)

# Each method against COSMIC, the neutral reference for "did it find a known
# signature".
cosmic_cmp <- rbind(
  cbind(method = "SignaturePPF", match_to_cosmic(fit$Signatures)),
  cbind(method = "TensorSignatures", match_to_cosmic(ts$signatures)))
write.csv(cosmic_cmp, file.path(DIR_TENSORSIG, "signature_cosmic_comparison.csv"),
          row.names = FALSE)
print(cosmic_cmp)

p_cos <- ggplot(cosmic_cmp, aes(method, cosine)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92") +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.8) +
  labs(x = NULL, y = "Best cosine similarity to COSMIC v3.4") +
  theme_bw()
ggsave(file.path(FIG_DIR, "TensorSignatures_cosine_to_cosmic.pdf"), p_cos,
       width = 4.5, height = 4)

################################################################################
# 7. Compare: chromatin-state effects
################################################################################
message("\n== 7. chromatin-state effects ==")
cmp <- compare_chromatin_effects(fit, TS_DIR, reference = REFERENCE_STATE)
write.csv(cmp, file.path(DIR_TENSORSIG, "chromatin_effect_comparison.csv"),
          row.names = FALSE)

ggsave(file.path(FIG_DIR, "TensorSignatures_chromatin_effects_pooled.pdf"),
       plot_chromatin_effects(cmp), width = 6, height = 5)
ggsave(file.path(FIG_DIR, "TensorSignatures_chromatin_effects_by_signature.pdf"),
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
# 8. Compare: regional mutation rate
#
#    The comparison that separates the two models. Aggregated into 1 Mb windows,
#    how close is each method's predicted burden to the observed one?
################################################################################
message("\n== 8. regional mutation rate ==")
rate <- compare_mutation_rate(dataChrom, fit, TS_DIR, window = 1e6)
write.csv(rate, file.path(DIR_TENSORSIG, "mutation_rate_windows.csv"),
          row.names = FALSE)

scores <- score_mutation_rate(rate)
write.csv(scores, file.path(DIR_TENSORSIG, "mutation_rate_scores.csv"),
          row.names = FALSE)
print(scores)

ggsave(file.path(FIG_DIR, "TensorSignatures_mutation_rate_along_genome.pdf"),
       plot_mutation_rate(rate), width = 12, height = 4)

message("\ndone: tables in ", DIR_TENSORSIG, "\n      figures in ", FIG_DIR)
