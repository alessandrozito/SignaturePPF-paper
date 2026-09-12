################################################################################
# Produces: Figure 3, and Figures S5, S8 and S12.
#           Also writes Figure S4, which FigureS4_signature_comparison.R now owns
#
# This file makes the plots for the denovo application.
#
# Usage:  Rscript R/Reproduce_figures_Application_denovo.R
#
# Everything is read from output/Application_denovo/, written by R/Application_denovo.R.
################################################################################

#---- Load packages
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(patchwork)

# Load the package
library(SignaturePPF)

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()
source(file.path(R_DIR, "Simulation_functions_main.R"))
source(file.path(R_DIR, "Simulation_functions.R"))

RERUN_COMPETITORS <- FALSE     # CompressiveNMF / SignatureAnalyzer, Figure S3

################################################################################
# Step 1 - Load the data and the output
################################################################################

data <- readRDS(PATH_ICGC2KB)
fit <- readRDS(file.path(DIR_DENOVO, "MCMCSolution.rds.gzip"))

# Relabel and add credible intervals
fit <- add_credible_intervals(fit, level = 0.95)
fit <- prune_signatures(fit)
fit <- relabel_by_mu(fit)

# Relabel to make coloring more pronounced in clustering
new_order <- c(1, 4, 2, 5, 3, 7, 6, 8, 9, 10)
fit <- relabel_by_mu(fit, order = new_order)

# Print the results
print(fit)

# Extract the posterior means and 95% quantiles
mu_order   <- fit$Mu
SigsMean   <- fit$Signatures
SigsLowCI  <- fit$lowCI$Signatures
SigsHighCI <- fit$highCI$Signatures
ThetaMean  <- fit$Thetas
BetasMean  <- fit$Betas
BetasLowCI <- fit$lowCI$Betas
BetasHighCI <- fit$highCI$Betas


################################################################################
# Step 2 - Figure 4
################################################################################

# Match to known COSMIC signatures
match_to_cosmic(SigsMean)

#---------------------------------------------------- Figure 4 panel a, b and c
# Plot the signatures, relevance weights, and regression coefficients
p_sigs  <- plot_Signatures(SigsMean, SigsLowCI, SigsHighCI)
p_mu    <- plot_Mu(fit, data = data) + theme(legend.position = "right")
p_betas <- plot_Betas(BetasMean, BetasLowCI, BetasHighCI)

p_sigs + p_mu + plot_spacer() + p_betas + plot_spacer()+
  plot_layout(widths = c(1.4, 0.2, 0.05, 1.5, 0.05))

ggsave(file.path(FIG_DIR, "Figure4_a_b_c_denovoPars.pdf"),
       width = 12.58, height = 5.82)

#---------------------------------------------------- Figure 4 panel d

# Cluster observations
clust <- cluster_samples(ThetaMean)
message("clusters: k = ", attr(clust, "k"), " (chosen by silhouette), sizes ",
        paste(table(clust), collapse = ", "))
print(round(attr(clust, "silhouette"), 4))

# Draw the plot
col_values <- c("darkblue", "#4959C7", "skyblue1", "lightblue1", "#F2AB67", "#F6C866",
                "#ED7470" , "brown", "#999999", "black")

plot_Theta(x = ThetaMean, clust = clust) +
  scale_fill_manual(values = col_values) +
  ylab(expression(paste("Total activities ", theta[kj]))) +
  theme(panel.grid.major.x = element_line(color = "lightgray"))

ggsave(file.path(FIG_DIR, "Figure4_d_denovoBaselines.pdf"),
       width = 12.34, height = 2.92)

# Matching patients and clusters with clinical data
data_clinical <- read_csv(file.path(OUTPUT_DIR, "PCAWG_clinical/clinical_BreastAdenoCA.csv"))
data_clust <- data_clinical %>%
  left_join(data.frame(clust) %>%
  rownames_to_column(var = "donor"),
  by = "donor")


# Compare clusters with cohorts
table(data_clust$clust, data_clust$project)
table(data_clust$clust, data_clust$grade)


################################################################################
# Step 3 - Figure 3 panel a: reconstruction at the Mb scale
################################################################################

LambdaPred <- reconstruct_lambda(fit, data$SignalTrack, data$CopyTrack)
LambdaPred <- rowSums(LambdaPred)
LambdaTrack <- data$gr_SignalTrack
LambdaTrack$LambdaPred <- LambdaPred

# Aggregate the 2 kb bins into 1 Mb windows
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6,
                      cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(LambdaTrack, hg19_Mb)
df_track <- as.data.frame(LambdaTrack) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred), .groups = "drop")

# Chromosome of each window
df_chrom <- as.data.frame(hg19_Mb) %>%
  mutate(region = seq_len(n()))

# Observed count per window, joined to the prediction and the chromosome
count_windows <- function(df_lambda) {
  as.data.frame(data$gr_Mutations) %>%
    mutate(region = subjectHits(findOverlaps(data$gr_Mutations, hg19_Mb))) %>%
    group_by(region) %>%
    summarise(n = n(), .groups = "drop") %>%
    left_join(df_lambda, by = "region") %>%
    left_join(dplyr::select(df_chrom, region, seqnames), by = "region") %>%
    mutate(
      chr_num = as.integer(gsub("chr", "", seqnames)),
      chrom_color_group = ifelse(is.na(chr_num %% 2), 1, chr_num %% 2),
      point_color_group = as.factor(chrom_color_group),
      line_color_group  = as.factor(chrom_color_group))
}

df_tumor_all <- count_windows(df_track)

# The alternating chromosome bands, shared by both track panels
chrom_bands <- as.data.frame(hg19_Mb) %>%
  mutate(region = seq_along(hg19_Mb)) %>%
  group_by(seqnames) %>%
  summarise(xmin = min(region), xmax = max(region) + 1, .groups = "drop") %>%
  arrange(seqnames) %>%
  mutate(fill_color = factor(seq_along(seqnames) %% 2))

plot_track <- function(df, linewidth = 0.7) {
  ggplot() +
    geom_rect(data = chrom_bands,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf,
                  fill = fill_color), alpha = 0.1) +
    scale_fill_manual(values = c("#4682B4", "antiquewhite")) +
    geom_point(data = df, size = 0.66,
               aes(x = region, y = n, color = point_color_group)) +
    scale_color_manual(values = c("#000D8B", "#93BAF1"), guide = "none") +
    ggnewscale::new_scale_color() +
    geom_line(data = df, linewidth = linewidth, alpha = 0.7,
              aes(x = region, y = Lambda, color = line_color_group,
                  group = seqnames)) +
    scale_color_manual(values = c("#CD2626", "tomato"), guide = "none") +
    theme_bw() +
    ylab("Number of mutations") +
    xlab("Genomic region (Mb)") +
    scale_x_continuous(expand = c(0.01, 0.02)) +
    scale_y_continuous(expand = c(0.03, 0.02)) +
    theme(legend.position = "none")
}

plot_reconstruction <- function(df, lim = c(0, 610)) {
  ggplot(df) +
    geom_point(aes(x = Lambda, y = n, color = point_color_group),
               alpha = 0.25, size = 0.5, shape = 20) +
    scale_color_manual(values = c("#000D8B", "#93BAF1"), guide = "none") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlim(lim) + ylim(lim) +
    geom_abline(slope = 1, intercept = 0, color = "#CD2626", linewidth = 0.7) +
    ylab("Observed mutations") +
    xlab("Predicted mutations")
}

p_Track <- plot_track(df_tumor_all)
p_reconstr <- plot_reconstruction(df_tumor_all)

################################################################################
# Step 4 - Figure 3 panel b: the same model without covariates
################################################################################

nocov_file <- file.path(DIR_DENOVO, "MAPSolution_noCovariates.rds.gzip")

if (file.exists(nocov_file)) {
  outNoCovs <- readRDS(nocov_file)
} else {
  message("\nfitting the covariate-free model...")
  outNoCovs <- SignaturePPF(
    data,
    sigs = NULL, sigs_fixed = FALSE, K = K_DENOVO,
    method = "map",
    prior = SignaturePPF_prior(),
    controls = SignaturePPF_control(maxiter = 500, tol = 1e-7,
                                    update_Betas = FALSE),
    init = SignaturePPF_init(
      Betas_start = matrix(0, ncol(data$SignalTrack), K_DENOVO)),
    seed = SEED, verbose = TRUE)
  saveRDS(outNoCovs, nocov_file, compress = "gzip")
}
outNoCovs <- relabel_by_mu(outNoCovs)
print(outNoCovs)

LambdaPredNoCovs <- rowSums(reconstruct_lambda(outNoCovs, data$SignalTrack,
                                               data$CopyTrack))
LambdaTrackNoCovs <- data$gr_SignalTrack
LambdaTrackNoCovs$LambdaPred <- LambdaPredNoCovs

df_trackNoCovs <- as.data.frame(LambdaTrackNoCovs) %>%
  mutate(region = subjectHits(over)) %>%
  group_by(region) %>%
  summarize(Lambda = sum(LambdaPred), .groups = "drop")

df_tumor_allNoCovs <- count_windows(df_trackNoCovs)

p_TrackNoCovs <- plot_track(df_tumor_allNoCovs, linewidth = 0.8)
p_reconstrNoCovs <- plot_reconstruction(df_tumor_allNoCovs)

(p_Track + p_reconstr) / (p_TrackNoCovs + p_reconstrNoCovs)
ggsave(file.path(FIG_DIR, "Figure3_recontructed_Mbscale.pdf"),
       width = 9.25, height = 4.65)

# The number the panel is making: how much of the regional structure each model
# recovers at the megabase scale.
regional <- data.frame(
  model = c("SignaturePPF", "copy number only"),
  rmse = c(sqrt(mean((df_tumor_all$Lambda - df_tumor_all$n)^2)),
           sqrt(mean((df_tumor_allNoCovs$Lambda - df_tumor_allNoCovs$n)^2))),
  cor = c(cor(df_tumor_all$Lambda, df_tumor_all$n),
          cor(df_tumor_allNoCovs$Lambda, df_tumor_allNoCovs$n)))
write.csv(regional, file.path(DIR_DENOVO, "regional_fit_Mb.csv"),
          row.names = FALSE)
print(regional)

################################################################################
# Step 5 - Figure S2: correlation across covariates
################################################################################

pdf(file.path(FIG_DIR, "FigureS2_covariate_correlation.pdf"),
    width = 7, height = 7)
# On the whole track
corrplot::corrplot(cor(data$SignalTrack), type = "lower",
                   addCoef.col = "gray25", number.cex = 0.7, diag = FALSE)
# At the mutations, which is the distribution the model actually sees
corrplot::corrplot(cor(mutation_covariates(data$gr_Mutations)), type = "lower",
                   addCoef.col = "gray25", number.cex = 0.7, diag = FALSE)
dev.off()

################################################################################
# Step 6 - Figure S3: PPF against CompressiveNMF and SignatureAnalyzer
################################################################################

mutMatrix <- getTotalMutations(data$gr_Mutations)
nmf_dir <- file.path(DIR_DENOVO, "BaselineNMF")
dir.create(nmf_dir, recursive = TRUE, showWarnings = FALSE)

f_compnmf <- file.path(nmf_dir, "out_CompNMF.rds")
f_sigan_kl <- file.path(nmf_dir, "out_SigAnalyzerL1KL.rds")
f_sigan_l2 <- file.path(nmf_dir, "out_SigAnalyzerL1WL2H.rds")

if (RERUN_COMPETITORS || !file.exists(f_compnmf)) {
  set.seed(10)
  BaselineNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix, K = 30,
                                                    alpha = 1.01, a = 1.01)
  saveRDS(BaselineNMF, f_compnmf)

  set.seed(10)
  out_SigAnalyzerL1KL <- sigminer::sig_auto_extract(nmf_matrix = t(mutMatrix),
                                                    method = "L1KL", cores = 20)
  saveRDS(out_SigAnalyzerL1KL, f_sigan_kl)

  set.seed(10)
  out_SigAnalyzerL1W.L2H <- sigminer::sig_auto_extract(nmf_matrix = t(mutMatrix),
                                                       method = "L1W.L2H",
                                                       cores = 20)
  saveRDS(out_SigAnalyzerL1W.L2H, f_sigan_l2)
}

BaselineNMF <- readRDS(f_compnmf)
out_SigAnalyzerL1KL <- readRDS(f_sigan_kl)
out_SigAnalyzerL1W.L2H <- readRDS(f_sigan_l2)

print(match_to_cosmic(BaselineNMF$Signatures))

match_to_cosmic(BaselineNMF$Signatures)
plot_Theta(BaselineNMF$Theta)

print(match_to_cosmic(out_SigAnalyzerL1KL$Signature.norm))

# Reconstruction of the count matrix under each model.
rmse <- c(
  SignaturePPF = sqrt(mean((fit$Signatures %*% fit$Thetas - mutMatrix)^2)),
  CompressiveNMF = sqrt(mean((BaselineNMF$Signatures %*% BaselineNMF$Theta -
                                mutMatrix)^2)),
  SignatureAnalyzer = sqrt(mean((out_SigAnalyzerL1KL$Signature.norm %*%
                                   out_SigAnalyzerL1KL$Exposure - mutMatrix)^2)))
print(round(rmse, 2))
write.csv(data.frame(model = names(rmse), rmse = as.numeric(rmse)),
          file.path(DIR_DENOVO, "count_matrix_rmse.csv"), row.names = FALSE)

# Matched one-to-one against the PPF signatures so the three panels line up
match_w_baseline <- match_MutSign(R_true = SigsMean,
                                  R_hat = BaselineNMF$Signatures)
colnames(match_w_baseline$R_hat) <- colnames(match_w_baseline$R_true)

match_w_sigAnalyzer <- match_MutSign(R_true = SigsMean,
                                     R_hat = out_SigAnalyzerL1KL$Signature.norm)
colnames(match_w_sigAnalyzer$R_hat) <- colnames(match_w_sigAnalyzer$R_true)

plot_Signatures(match_w_baseline$R_true) +
  plot_Signatures(match_w_baseline$R_hat) +
  plot_Signatures(match_w_sigAnalyzer$R_hat)
ggsave(file.path(FIG_DIR, "Breast_suppl_Signatures_comparison.pdf"),
       width = 13.42, height = 5.64)

################################################################################
# Step 7 - Figure S4: effective sample sizes and the log posterior trace
################################################################################

keep <- kept_draw_index(fit)
all_draws <- seq_along(keep)

sigs <- colnames(fit$Signatures)

ThetaEffects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$THETAchain[keep, sigs, , drop = FALSE], all_draws)
REffects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$SIGSchain[keep, , sigs, drop = FALSE], all_draws)
BetasEffects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$BETASchain[keep, , sigs, drop = FALSE], all_draws)
MuEffects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$MUchain[keep, sigs, drop = FALSE], all_draws)

effects_df <- tibble(
  value = c(REffects, ThetaEffects, BetasEffects, MuEffects),
  quantity = factor(rep(c("R", "Theta", "B", "mu"),
                        times = c(length(REffects), length(ThetaEffects),
                                  length(BetasEffects), length(MuEffects))),
                    levels = c("R", "Theta", "B", "mu")))

p_effective <- ggplot(effects_df, aes(x = "", y = value)) +
  geom_boxplot() +
  facet_wrap(~ quantity, nrow = 1, scales = "free") +
  labs(x = "", y = "Effective Sample Size") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# logpost_every = 10 leaves the other nine iterations NA, so the trace is
# subsampled rather than complete
lp <- fit$MCMCchain$logPostchain[keep]
lp_kept <- which(!is.na(lp))
p_logPost <- ggplot(data.frame(iter = lp_kept * fit$controls$thin,
                               logPosterior = lp[lp_kept])) +
  geom_line(aes(x = iter, y = logPosterior)) +
  theme_bw() +
  facet_wrap(~ "logposterior trace", nrow = 1) +
  xlab("post-burnin iteration")

p_ess <- p_effective + p_logPost + plot_layout(nrow = 1, widths = c(2, 1))
ggsave(file.path(FIG_DIR, "FigureS4_ess_denovo.pdf"), p_ess,
       width = 11.68, height = 3.5)

################################################################################
# Step 8 - Sensitivity of the MAP solution to Kmax and to the beta prior
#
################################################################################
message("\n== sensitivity of the MAP solution ==")

# Signatures in the MAP with K = 12
ref_map <- relabel_by_mu(prune_signatures(
  readRDS(file.path(DIR_DENOVO, "MAPSolution.rds.gzip")), verbose = FALSE))

mm_ref <- match_MutSign(R_true = fit$Signatures, R_hat = ref_map$Signatures)
ref_map <- relabel_by_mu(ref_map,
                         order = mm_ref$match[seq_len(ncol(fit$Signatures))])
message("MAP vs reported solution, cosine per signature:")
print(round(stats::setNames(
  colSums(fit$Signatures * ref_map$Signatures) /
    sqrt(colSums(fit$Signatures^2) * colSums(ref_map$Signatures^2)),
  colnames(fit$Signatures)), 3))

# All scenarios
scenarios <- c(K15 = "MAP_K15.rds.gzip",
               K20 = "MAP_K20.rds.gzip",
               K12_c0_10 = "MAP_K12_c0_10.rds.gzip",
               K12_c0_1 = "MAP_K12_c0_1.rds.gzip",
               K20_c0_10 = "MAP_K20_c0_10.rds.gzip",
               K20_c0_1 = "MAP_K20_c0_1.rds.gzip")
absent <- !file.exists(file.path(DIR_SENSITIVITY, scenarios))
if (any(absent)) {
  message("not yet fitted, left out: ",
          paste(names(scenarios)[absent], collapse = ", "))
}
scenarios <- scenarios[!absent]

# A reference signature counts as recovered at this cosine or above.
SENS_COS_CUTOFF <- 0.9

sens <- do.call(rbind, lapply(names(scenarios), function(nm) {
  raw <- readRDS(file.path(DIR_SENSITIVITY, scenarios[[nm]]))
  fit <- relabel_by_mu(prune_signatures(raw, verbose = FALSE))

  # Matched one-to-one against the reference, then scored on the pairs that are
  # real: match_MutSign() pads the smaller set to make the assignment square and
  # zeroes the padding on the way out, and a zero column's cosine is NaN.
  mm <- match_MutSign(R_true = ref_map$Signatures, R_hat = fit$Signatures)
  paired <- colSums(mm$R_true) > 0 & colSums(mm$R_hat) > 0
  cs <- colSums(mm$R_true[, paired, drop = FALSE] * mm$R_hat[, paired, drop = FALSE]) /
    sqrt(colSums(mm$R_true[, paired, drop = FALSE]^2) *
           colSums(mm$R_hat[, paired, drop = FALSE]^2))

  data.frame(
    scenario = nm,
    # `K_fitted` is set by prune_signatures(), so it lives on the pruned copy -
    # these fits were run with prune_solution = FALSE and `raw` has no such field.
    Kmax = fit$K_fitted,
    c0 = raw$prior$c0, d0 = raw$prior$d0,
    n_selected = fit$K,
    min_mu = min(fit$Mu),
    matched = sum(cs >= SENS_COS_CUTOFF),
    mean_cosine = mean(cs),
    iterations = as.integer(raw$MAPsolution$iter),
    sec_per_iter = as.numeric(raw$runtime, units = "secs") /
      as.integer(raw$MAPsolution$iter),
    hours = as.numeric(raw$runtime, units = "hours"),
    row.names = NULL, stringsAsFactors = FALSE)
}))

reference_row <- data.frame(
  scenario = "reference", Kmax = ref_map$K_fitted,
  c0 = ref_map$prior$c0, d0 = ref_map$prior$d0,
  n_selected = ref_map$K,
  min_mu = min(ref_map$Mu),
  matched = ref_map$K, mean_cosine = 1,   # matched against itself
  iterations = as.integer(ref_map$MAPsolution$iter),
  sec_per_iter = as.numeric(ref_map$runtime, units = "secs") /
    as.integer(ref_map$MAPsolution$iter),
  hours = as.numeric(ref_map$runtime, units = "hours"),
  row.names = NULL, stringsAsFactors = FALSE)

sens <- rbind(reference_row, sens)
write.csv(sens, file.path(DIR_SENSITIVITY, "sensitivity_summary.csv"),
          row.names = FALSE)
print(format(sens, digits = 3), row.names = FALSE)


stopifnot(all(sens$matched > 0))

converged <- vapply(c(list(ref_map), lapply(scenarios, function(f)
  readRDS(file.path(DIR_SENSITIVITY, f)))),
  function(x) x$MAPsolution$iter < x$controls$maxiter, logical(1))
if (!all(converged)) {
  warning("scenario(s) stopped at maxiter rather than converging: ",
          paste(sens$scenario[!converged], collapse = ", "), call. = FALSE)
}

################################################################################
# Step 8b - Every scenario against the reference, signature by signature
#
################################################################################

live <- function(f) relabel_by_mu(prune_signatures(f, verbose = FALSE))
fits <- c(list(reference = ref_map),
          lapply(scenarios, function(f) live(readRDS(file.path(DIR_SENSITIVITY, f)))))

# ---- the union of signatures, reference first
anchor <- ref_map$Signatures
for (nm in names(scenarios)) {
  S <- fits[[nm]]$Signatures
  mm <- match_MutSign(R_true = anchor, R_hat = S)
  taken <- mm$match[seq_len(ncol(anchor))]
  new <- setdiff(seq_len(ncol(S)), taken[taken <= ncol(S)])
  if (length(new)) {
    add <- S[, new, drop = FALSE]
    colnames(add) <- sprintf("Extra%02d",
                             ncol(anchor) - ncol(ref_map$Signatures) + seq_along(new))
    anchor <- cbind(anchor, add)
  }
}
n_ref <- ncol(ref_map$Signatures)
message("union: ", ncol(anchor), " rows = ", n_ref, " reference + ",
        ncol(anchor) - n_ref, " found only by a scenario")

# ---- each scenario placed onto those rows
onto_anchor <- function(f) {
  mm <- match_MutSign(R_true = anchor, R_hat = f$Signatures)
  partner <- mm$match[seq_len(ncol(anchor))]
  partner[partner > ncol(f$Signatures)] <- NA_integer_   # a padding column

  cosine <- rep(NA_real_, ncol(anchor))
  ok <- which(!is.na(partner))
  cosine[ok] <- vapply(ok, function(i) {
    a <- anchor[, i]; b <- f$Signatures[, partner[i]]
    sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  }, numeric(1))

  assigned <- df_assign(f, data)
  S <- matrix(0, nrow(anchor), ncol(anchor), dimnames = dimnames(anchor))
  S[, ok] <- f$Signatures[, partner[ok], drop = FALSE]

  list(
    table = data.frame(
      signature = colnames(anchor),
      partner = colnames(f$Signatures)[partner],
      cosine = cosine,
      present = !is.na(partner),
      mu = as.numeric(f$Mu)[partner],
      m = assigned$m[match(colnames(f$Signatures)[partner], assigned$best_sig)],
      row.names = NULL, stringsAsFactors = FALSE),
    Signatures = S)
}

matches <- lapply(fits, onto_anchor)

match_tbl <- do.call(rbind, lapply(names(matches), function(nm)
  cbind(scenario = nm, matches[[nm]]$table)))
write.csv(match_tbl, file.path(DIR_SENSITIVITY, "signature_matching.csv"),
          row.names = FALSE)
print(format(match_tbl, digits = 3), row.names = FALSE)

#---------------------------------------------- panel labels for the scenarios
scenario_label <- sprintf("(%g, %g, %g)", sens$Kmax, sens$c0, sens$d0)
names(scenario_label) <- sens$scenario
scenario_label <- scenario_label[names(matches)]   # panel order
stopifnot(!anyNA(scenario_label))

#---------------------------------------------- relevance weights, one grid
mu_grid <- match_tbl[match_tbl$present, ]
mu_grid$scenario <- factor(scenario_label[mu_grid$scenario],
                           levels = unname(scenario_label))
mu_grid$signature <- factor(mu_grid$signature, levels = colnames(anchor))

p_mu_sens <- ggplot(mu_grid, aes(x = scenario, y = 1, size = mu, fill = m)) +
  geom_point(colour = "black", stroke = 0.7, shape = 21) +
  facet_wrap(~ signature, ncol = 1, strip.position = "left") +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_size(name = expression(mu[k]), range = c(3, 12)) +
  scale_fill_gradientn(
    name = "N. mutations",
    colours = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0,
                                   size = 10, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.position = "right",
        panel.spacing = grid::unit(0.03, "lines"))

ggsave(file.path(FIG_DIR, "Sensitivity_mu_all_scenarios.pdf"), p_mu_sens,
       width = 2.3 + 0.60 * length(matches), height = 7.10)

#------------------------------------- how COSMIC-like is each row on average?
cosmic_rows <- do.call(rbind, lapply(names(matches), function(nm) {
  S <- matches[[nm]]$Signatures
  # onto_anchor() leaves an all-zero column wherever a scenario has no partner
  # for that row. Those are absences, not signatures, and must not be scored.
  keep <- which(colSums(S) > 0)
  if (!length(keep)) return(NULL)
  m <- match_to_cosmic(S[, keep, drop = FALSE])
  data.frame(scenario = nm, signature = colnames(S)[keep],
             best_match = m$best_match, cosine = m$cosine,
             row.names = NULL, stringsAsFactors = FALSE)
}))

cosmic_tbl <- do.call(rbind, lapply(colnames(anchor), function(sg) {
  d <- cosmic_rows[cosmic_rows$signature == sg, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  tb <- sort(table(d$best_match), decreasing = TRUE)
  data.frame(signature = sg,
             n_scenarios = nrow(d),
             mean_cosine = mean(d$cosine),
             sd_cosine = if (nrow(d) > 1) stats::sd(d$cosine) else 0,
             # the COSMIC signature matched most often across the scenarios;
             # n_match records how often, so a row that is not consistently
             # assigned to one COSMIC signature is visible as such.
             best_match = names(tb)[1],
             n_match = as.integer(tb[1]),
             row.names = NULL, stringsAsFactors = FALSE)
}))
cosmic_tbl$signature <- factor(cosmic_tbl$signature, levels = colnames(anchor))
write.csv(cosmic_tbl, file.path(DIR_SENSITIVITY, "signature_cosmic_cosine.csv"),
          row.names = FALSE)
print(format(cosmic_tbl, digits = 3), row.names = FALSE)

p_cosmic <- ggplot(cosmic_tbl, aes(x = mean_cosine, y = 1)) +
  geom_col(orientation = "y", fill = "grey80", width = 0.5) +
  geom_text(aes(label = sprintf("%.3f", mean_cosine)),
            x = 0.04, hjust = 0, size = 2.9, colour = "black") +
  geom_text(aes(label = best_match),
            x = 0.96, hjust = 1, size = 2.9, colour = "grey30") +
  facet_grid(signature ~ ., switch = "y") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  #ggtitle("Cosine to\nclosest COSMIC") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10, hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text.y.left = element_blank(),
        panel.spacing.y = grid::unit(0, "lines"))

#---------------------------------------------- the spectra, side by side
spectra_panels <- lapply(names(matches), function(nm) {
  p <- plot_Signatures(matches[[nm]]$Signatures) +
    ggtitle(scenario_label[[nm]]) +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          axis.text.x = element_blank(),
          strip.text.x = element_blank())
  # The signature labels only need to appear once, down the left-hand edge.
  if (nm != names(matches)[1]) p <- p + theme(strip.text.y.left = element_blank())
  p
})

# The cosine column is the last panel, narrow, and shares the rows.
p_spectra <- Reduce(`+`, spectra_panels) + p_cosmic +
  plot_layout(nrow = 1, widths = c(rep(1, length(spectra_panels)), 0.80))

#p_mu_sens + p_spectra + plot_layout(widths = c(0.25, 1))

ggsave(file.path(FIG_DIR, "Sensitivity_signatures_comparison.pdf"), p_spectra,
       width = 16.02, height = 5.43, limitsize = FALSE)

################################################################################
# Step 9 - Additional results
################################################################################

#---- Summary of the MAP restarts: cost, and whether they agree
map_starts <- read.csv(file.path(DIR_DENOVO, "MapSolutions", "map_starts.csv"))
map_starts$seconds_per_iteration <- map_starts$minutes * 60 / map_starts$iterations
print(map_starts)

maps <- lapply(seq_len(nrow(map_starts)), function(i) {
  relabel_by_mu(prune_signatures(
    readRDS(file.path(DIR_DENOVO, "MapSolutions",
                      sprintf("MAPSolution_start%02d.rds.gzip", i))),
    verbose = FALSE))
})

# Do the restarts find the same signatures? Matched pairwise, then the mean
# cosine on the diagonal - one number per pair.
cosine_diag <- function(a, b) {
  m <- match_MutSign(R_true = a, R_hat = b)
  paired <- colSums(m$R_true) > 0 & colSums(m$R_hat) > 0
  cs <- colSums(m$R_true[, paired, drop = FALSE] * m$R_hat[, paired, drop = FALSE]) /
    sqrt(colSums(m$R_true[, paired, drop = FALSE]^2) *
           colSums(m$R_hat[, paired, drop = FALSE]^2))
  c(matched = sum(paired), mean_cosine = mean(cs), min_cosine = min(cs))
}
pairs <- utils::combn(length(maps), 2)
scores <- t(apply(pairs, 2, function(p)
  cosine_diag(maps[[p[1]]]$Signatures, maps[[p[2]]]$Signatures)))
agreement <- data.frame(
  pair = paste0(pairs[1, ], "-", pairs[2, ]),
  n_signatures = apply(pairs, 2, function(p)
    paste0(maps[[p[1]]]$K, "/", maps[[p[2]]]$K)),
  scores)
print(agreement)
write.csv(agreement, file.path(DIR_DENOVO, "map_restart_agreement.csv"),
          row.names = FALSE)

#---- Introduction: chromosome 1p against 1q
chr1_arms_gr <- GRanges(
  seqnames = c("chr1", "chr1"),
  ranges = IRanges(start = c(1, 124535434),
                   end = c(121535433, 249250621)),
  arm = c("1p", "1q"))

df_chrom1pq <- as.data.frame(data$gr_Mutations)
df_chrom1pq$region <- NA
over1pq <- findOverlaps(data$gr_Mutations, chr1_arms_gr)
df_chrom1pq$region[queryHits(over1pq)] <- subjectHits(over1pq)

df_copy1pq <- as.data.frame(data$gr_CopyTrack)
df_copy1pq$region <- NA
overCopy1pq <- findOverlaps(data$gr_CopyTrack, chr1_arms_gr)
df_copy1pq$region[queryHits(overCopy1pq)] <- subjectHits(overCopy1pq)

arm_summary <- df_chrom1pq %>%
  filter(!is.na(region)) %>%
  group_by(region) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(df_copy1pq %>%
              filter(!is.na(region)) %>%
              dplyr::select(region, all_of(colnames(data$CopyTrack))) %>%
              drop_na() %>%
              gather(key = "key", value = "value", -region) %>%
              group_by(region) %>%
              summarise(Copies = 2 * mean(value), .groups = "drop"),
            by = "region") %>%
  mutate(arm = chr1_arms_gr$arm[region])
print(arm_summary)

message("\ndone: figures in ", FIG_DIR)
