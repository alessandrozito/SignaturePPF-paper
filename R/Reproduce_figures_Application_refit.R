################################################################################
# Produces: Figure 4, and Figure S6
#
# This pre-specified signatures application in the paper.
################################################################################

#---- Load packages
library(tidyverse)
library(BSgenome)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(foreach)
library(patchwork)
library(doParallel)
library(GenomicRanges)
library(AnnotationDbi)

# Load the package
library(SignaturePPF)

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()
source(file.path(R_DIR, "Simulation_functions_main.R"))

######################################################################
# Plotting functions

#--- Useful functions
## mutation_covariates() now lives in Utils_functions.R: the de novo figures need
## it too, and two copies is one copy too many.

Compute_mutation_Probs <- function(gr_Mutations, R, Theta, Betas){
  X <- mutation_covariates(gr_Mutations)
  Probs <- R[gr_Mutations$channel, ] * t(Theta[, gr_Mutations$sample]) * exp(X %*% Betas)
  Probs <- t(apply(Probs, 1, function(x) x/sum(x)))
  colnames(Probs) <- colnames(Betas)
  return(Probs)
}

get_signal_regions <- function(chr, region_start, region_end, gr_SignalTrack) {
  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  return(gr_subset)
}

find_mutations_in_region <- function(data,
                                     i, j,
                                     chr,
                                     region_start,
                                     region_end,
                                     R,
                                     Theta,
                                     Betas){
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  gr_mut_filter <- data$gr_Mutations[data$gr_Mutations$sample == j & data$gr_Mutations$channel == i]
  overlaps <- findOverlaps(gr_region, gr_mut_filter)
  positions <- start(gr_mut_filter[subjectHits(overlaps)])

  # Now, find the value of the signal at the mutations
  X <- as.matrix(mcols(get_signal_regions(chr, positions, positions + 1, data$gr_SignalTrack)))[, -1]
  intensity <- t(R[i, ] * Theta[, j])[rep(1, length(positions)), ] * exp(X %*% Betas)
  max_int <- apply(intensity, 1, which.max)
  df_pos <- data.frame("positions" = positions, "Sig" = colnames(Betas)[max_int])
  return(df_pos)
}

calculate_ChannelProbs_region <- function(chr,
                                          region_start, region_end,
                                          R,
                                          Theta,
                                          Betas,
                                          j = 1,
                                          gr_SignalTrack,
                                          CopyTrack) {


  # Create the region
  gr_region <- GRanges(chr, IRanges(region_start + 1, region_end))
  # Subset the signalTrack into relevant regions
  overlaps <- findOverlaps(gr_region, gr_SignalTrack)
  gr_subset <- gr_SignalTrack[subjectHits(overlaps)]
  # Estimate the loadings in each subset
  SignalTrack_subset <- as.matrix(mcols(gr_subset)[, -1])
  CopyTrack_subset <- CopyTrack[subjectHits(overlaps), j]
  ExpSolBetas <- exp(SignalTrack_subset %*% Betas)
  # Big Storage
  ChannelProbs <- array(NA, dim = c(nrow(ExpSolBetas), ncol(ExpSolBetas), 96),
                        dimnames = list(NULL, colnames(ExpSolBetas), rownames(R)))

  ExpSolBetas_copy <- apply(ExpSolBetas, 2, function(x) x * CopyTrack_subset)
  for(i in 1:96){
    SigTheta <- t(R[i, ] * Theta[, j])
    bigProd <- SigTheta[rep(1, nrow(ExpSolBetas_copy)), ] * ExpSolBetas_copy
    sig_probs <-  bigProd#/rowSums(bigProd)
    ChannelProbs[, , i] <- sig_probs
  }
  ChannelProbs
}

################################################################################
# Step 1 - Load the data and the output
################################################################################

# Load the data
data <- readRDS(PATH_ICGC2KB)

fit <- readRDS(file.path(DIR_REFIT, "MCMCSolution.rds.gzip"))
results <- readRDS(file.path(DIR_REFIT, "resultsMCMC_refit.rds.gzip"))

#-- Signatures
Sigs_to_use <- SIGS_TO_USE
CosmicSigs <- SignaturePPF::COSMIC_v3.4_SBS96_GRCh37[, Sigs_to_use]

#----------------------------------- Figure 5 panel a
# Regression coefficients
p_Beats_sigs_fixed <- plot_Betas(results$Betas$mean, results$Betas$lowCI, results$Betas$highCI)

# Relevance weights
p_mu_sigs_fixed <- plot_Mu(fit, data = data) + theme(legend.position = "right")

#----------------------------------- Figure 5 panel b

#--- Run the baseline CompressiveNMF
mutMatrix <- getTotalMutations(data$gr_Mutations)
BaseNMF <- CompressiveNMF::CompressiveNMF_map(mutMatrix,
                                              K = 0,
                                              S = 1e6 * CosmicSigs + 1,
                                              alpha = 1.01, a = 1.01)
dimnames(BaseNMF$mapOutput$R) <- dimnames(CosmicSigs)
colnames(BaseNMF$mapOutput$Theta) <- colnames(mutMatrix)
rownames(BaseNMF$mapOutput$Theta) <- colnames(CosmicSigs)


# Probabilities for each mutation in PPF
MutProbs <- Compute_mutation_Probs(data$gr_Mutations,
                                   R = fit$Signatures,
                                   Theta = fit$Baseline,
                                   Betas = results$Betas$mean)

# Probabilities for each mutation in NMF
Betas_zero <- matrix(0,
                     nrow = nrow(results$Betas$mean),
                     ncol = ncol(results$Betas$mean),
                     dimnames = dimnames(results$Betas$mean))
MutProbs_fixed <- Compute_mutation_Probs(data$gr_Mutations,
                                         R = BaseNMF$mapOutput$R,
                                         Theta = BaseNMF$mapOutput$Theta,
                                         Betas = Betas_zero)

# Signature with the highest probability in both cases
best_sigs_ppf <- Sigs_to_use[apply(MutProbs, 1, which.max)]
best_sigs_ppf <- factor(best_sigs_ppf, levels = Sigs_to_use)

best_sigs_nmf <- Sigs_to_use[apply(MutProbs_fixed, 1, which.max)]
best_sigs_nmf <- factor(best_sigs_nmf, levels = Sigs_to_use)

df_probs <- data.frame("Best_ppf" = best_sigs_ppf,
                       "Prob_ppf" = apply(MutProbs, 1, max),
                       "Best_nmf" = best_sigs_nmf,
                       "Prob_nmf" = apply(MutProbs_fixed, 1, max))

# Plot the confusion matrix between PPF and baseline NMF
p_confusion <- df_probs %>%
  group_by(Best_ppf, Best_nmf) %>%
  summarise(
    count = n(),
    mean_prob_ppf = mean(Prob_ppf),
    .groups = "drop"
  ) %>%
  complete(Best_ppf, Best_nmf) %>%
  mutate(AltCol = ifelse((as.numeric(Best_ppf) %% 2) == 0, "gray93", "gray97")) %>%
  ggplot(aes(x = Best_ppf, y = Best_nmf)) +
  geom_tile(aes(fill = AltCol), color = "white", width = 0.95, height = 0.95) +
  scale_fill_manual(values = c("gray93", "gray97"), guide = "none") +
  geom_point(aes(size = count, color = mean_prob_ppf), na.rm = TRUE) +
  scale_size_binned(
    range = c(0.5, 10),
    breaks = c(10, 100, 1000, 5000, 10000, 40000, 100000, 140000),
    name = "Count"
  ) +
  ggplot2::scale_color_gradientn(
    name = "Average\nprobability",
    colours = rev(c("#fde725", "#b5de2b", "#6ece58", "#35b779", "#1f9e89",
                    "#26828e", "#31688e", "#3e4989", "#482878", "#440154"))) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(df_probs$Best_nmf))) +
  labs(x = "best_sigs_ppf", y = "best_sigs_nmf", size = "Count") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0),
    panel.grid = element_blank(),
    aspect.ratio = 1,
    legend.position = "right",        # <--- move legend top
    legend.box = "horizontal"       # <--- stack legends horizontally
  )


#----- Save the outputs
ggsave(plot = p_mu_sigs_fixed, filename = file.path(FIG_DIR, "Figure5_a_mu.pdf"),
       width = 8.42, height = 5.56)

ggsave(plot = (p_Beats_sigs_fixed + theme(aspect.ratio = 1)) + p_confusion,
       filename = file.path(FIG_DIR, "Figure5_a_b.pdf"),
       width = 11.68, height = 5.05)


#----------------------------------- Figure 5 panel c
i <- "A[C>A]C"
j <- "DO220823"
chr <- "chr1"
region_start <- 150e6
region_end   <- 160e6

ChannelProbs <- calculate_ChannelProbs_region(chr,
                                              region_start,
                                              region_end,
                                              fit$Signatures,
                                              fit$Baseline,
                                              results$Betas$mean,
                                              j = j,
                                              data$gr_SignalTrack,
                                              data$CopyTrack)

pr_const <- 2000 *BaseNMF$mapOutput$R[i, ] * BaseNMF$mapOutput$Theta[, j]/(sum(data$gr_SignalTrack$bin_weight))
df_const <- data.frame(t(pr_const)[rep(1,nrow(ChannelProbs)), ]) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1)) * 2000 + region_start,
         method = "Standard NMF")

df_mod <- data.frame(ChannelProbs[, , i])
colnames(df_mod) <- colnames(BaseNMF$mapOutput$R)

df_pos <- find_mutations_in_region(data, i, j, chr,
                                   region_start, region_end,
                                   fit$Signatures,
                                   fit$Baseline,
                                   results$Betas$mean)

#------- Top panel
window <- 10
options(scipen=99)
p_intensity <- df_mod %>%
  mutate(across(where(is.numeric), ~ zoo::rollmean(.x, fill = NA, k = window,
                                                   align = "center"))) %>%
  mutate(region = (0:(nrow(df_mod) - 1))*2000 + region_start,
         method = "PoissonProcess") %>%
  bind_rows(df_const) %>%
  gather(key = "Sig", value = "Prob", -region, -method) %>%
  filter(Sig %in% c("SBS3", "SBS40a", "SBS8")) %>%
  ggplot() +
  geom_line(aes(x = region, y = Prob, color = Sig, linetype = method, alpha = method),
            linewidth = 0.6) +
  scale_color_manual(name = "Signature",
                     values = c("red", "antiquewhite3", "blue")) +
  scale_alpha_manual(values = c(1, 1))+
  ylab("Expected mutations in 2kb region") +
  facet_wrap(~paste0("Mutation type ", i, " - Patient ", j)) +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  geom_rug(data = df_pos,
           aes(x = positions, color = Sig),
           sides = "b", # "b" for bottom
           inherit.aes = FALSE,
           length = grid::unit(0.035, "npc"),
           linewidth = 1.2)+
  theme(axis.title.x = element_blank())

#------- Covariates
subset_tracks <- get_signal_regions(chr, region_start, region_end, data$gr_SignalTrack)
subset_tracks$bin_weight <- NULL
df_tracks <- as.data.frame(mcols(subset_tracks)) %>%
  mutate(region = (0:(nrow(ChannelProbs) - 1))*2000 + region_start) %>%
  gather(key = "covariate", value = "score", -region)

p_track <- ggplot(df_tracks) +
  geom_raster(aes(x = region + 1000, y = covariate, fill = score)) +
  #facet_wrap(.~ pos2, nrow = 1) +
  scale_fill_gradient2(name = "Signal track\n(standardized)",
                       mid = "white", low = "#2166AC", high = "#B2182B", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish)+
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  #theme(axis.title = element_blank())
  ylab("Genomic covariate")+
  xlab(paste0("Genomic regions in ",  chr, " (2kb)"))
p_panels <- p_intensity/p_track + plot_layout(heights = c(1.6, 1))
ggsave(plot = p_panels, filename = file.path(FIG_DIR, "Figure5_c_.png"),
       width = 11.68, height = 5.91)


################################################ Figure S5
# Indexes to keep - postburnin
keep <- kept_draw_index(fit)

# ---- Recalculate mu for this solution (post-burnin)
mu <- colMeans(fit$MCMCchain$MUchain[keep, , drop = FALSE])
# ---- Extract the full chain of Thetas (total activities)
ThetaChain_adj <- fit$MCMCchain$THETAchain[keep, mu > 0.01, , drop = FALSE]

# ---- Calculate effective sample sizes
all_draws <- seq_len(dim(ThetaChain_adj)[1])
ThetaEffects <- get_PosteriorEffectiveSize(ThetaChain_adj, all_draws)
BetasEffects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$BETASchain[keep, , mu > 0.01, drop = FALSE], all_draws)
Sigma2Effects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$SIGMA2chain[keep, mu > 0.01, drop = FALSE], all_draws)
MuEffects <- get_PosteriorEffectiveSize(
  fit$MCMCchain$MUchain[keep, mu > 0.01, drop = FALSE], all_draws)

# ---- Make tidy data frame
effects_df <- tibble(
  value = c(ThetaEffects, BetasEffects, MuEffects),
  quantity = factor(rep(
    c("Theta([0, T))", "B", "mu"),
    times = c(length(ThetaEffects), length(BetasEffects), length(MuEffects))
  ),
  levels = c("Theta([0, T))", "B", "mu"))
)

# ---- Effective sample size boxplot
p_effective <- ggplot(effects_df, aes(x = "", y = value)) +
  geom_boxplot() +
  #geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  facet_wrap(~ quantity, nrow = 1, scales = "free") +
  labs(x = "", y = "Effective Sample Size") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ---- Log posterior trace plot
# logpost_every = 10 leaves the other nine sweeps NA, so the trace is subsampled
# rather than complete; the x axis is the sweep it was evaluated at.
lp <- fit$MCMCchain$logPostchain[keep]
lp_kept <- which(!is.na(lp))
p_logPost <- ggplot(
  data = data.frame(
    iter = lp_kept * fit$controls$thin,
    logPosterior = lp[lp_kept]
  )
) +
  geom_line(aes(x = iter, y = logPosterior)) +
  theme_bw() +
  facet_wrap(~ "logposterior trace", nrow = 1) +
  xlab("post-burnin iteration")

p_ess_fixed <- p_effective + p_logPost + plot_layout(nrow = 1, widths = c(2, 1))

# The original only displayed this one; saved here so a batch run keeps it.
ggsave(plot = p_ess_fixed, filename = file.path(FIG_DIR, "FigureS5_ess_refit.pdf"),
       width = 11.68, height = 3.5)

message("done: figures in ", FIG_DIR)
