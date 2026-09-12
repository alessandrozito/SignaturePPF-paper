################################################################################
# Produces: Figure S10
#
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

invisible(SignaturePPF_validate(data80))
invisible(SignaturePPF_validate(dataICGC))

stopifnot(identical(colnames(data80$SignalTrack), colnames(dataICGC$SignalTrack)),
          nrow(data80$SignalTrack) == nrow(dataICGC$SignalTrack))

CosmicSigs <- COSMIC_v3.4_SBS96_GRCh37[, SIGS_TO_USE]

# Describe the number of mutations and overlaps
sort(table(data80$gr_Mutations$sample))
quantile(table(data80$gr_Mutations$sample))
hist(table(data80$gr_Mutations$sample), breaks = 30)
sum(table(data80$gr_Mutations$sample) > 9000)

################################################################################
# 2. Visualize the mutations in both datasets
################################################################################
#------------------------------------------------------------------
# Plot the mutation distribution along the genome
genome <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]
hg19_Mb <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)

#---- Count the number of mutations in the dataset per megabase
df_tumor80 <- as.data.frame(data80$gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(data80$gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  mutate(data = "Breast80")

df_tumorICGC <- as.data.frame(dataICGC$gr_Mutations) %>%
  mutate(region = subjectHits(findOverlaps(dataICGC$gr_Mutations, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n()) %>%
  mutate(data = "BreastICGC")

df_tumor_all <- rbind(df_tumor80, df_tumorICGC)

# Plot the Mb mutations side by side
p_breast_both <- ggplot() +
  theme_bw(base_size = 11) +
  geom_rect(data = as.data.frame(hg19_Mb) %>%
              mutate(region = 1:length(hg19_Mb)) %>%
              group_by(seqnames) %>%
              summarise(xmin = min(region), xmax = max(region) + 1) %>%
              arrange(seqnames) %>%
              mutate(fill_color = factor(seq_along(seqnames) %% 2)),
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
            alpha = 0.2) +
  scale_fill_manual(values = c("#4682B4", "antiquewhite2" )) +
  geom_point(data = df_tumor_all,
             aes(x = region, y = n), size = 1, shape = 21, color = "gray25") +
  xlab("Genomic region (1Mb)") +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  scale_x_continuous(expand = c(0.01, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  facet_wrap(.~data, nrow = 2, scales ="free_y")

# The number of mutations in tumor80 dataset is almost half of ICGC
df_chroms <- as.data.frame(hg19_Mb) %>%
  mutate(region = 1:length(hg19_Mb)) %>%
  group_by(seqnames) %>%
  summarise(xmin = min(region), xmax = max(region) + 1) %>%
  arrange(seqnames) %>%
  mutate(fill_color = factor(seq_along(seqnames) %% 2)) %>%
  rowwise() %>%
  mutate(region = list(xmin:(xmax - 1))) %>%
  unnest(region) %>%
  ungroup()

df_tumor_all %>%
  left_join(df_chroms, by = "region") %>%
  dplyr::select(-xmin, -xmax) %>%
  pivot_wider(id_cols = c(region, seqnames, fill_color),
              names_from = data,
              values_from = n) %>%
  ggplot(aes(x = Breast80, y = BreastICGC)) +
  geom_point(aes(color = fill_color), alpha = 0.5, size = 0.5) +
  scale_color_manual(values = c("#4682B4", "darkorange"))+
  theme_bw() +
  xlab("Mb mutations in Breast80") +
  ylab("Mb mutations in BreastICGC")+
  theme(legend.position = "none")

mutMatrix80 <- getTotalMutations(data80$gr_Mutations)
mutMatrixICGC <- getTotalMutations(dataICGC$gr_Mutations)

p_total_muts <- CompressiveNMF::plot_SBS_signature(cbind("BreastICGC" = rowSums(mutMatrixICGC), "Breast80" = rowSums(mutMatrix80))) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 8),
        panel.grid.major.y  = element_line(linewidth = 0.1, color = "gray"))

p_all <- p_breast_both + p_total_muts + plot_layout(widths = c(2,1))
ggsave(file.path(FIG_DIR, "Breast80_BreastICGC_descripion.pdf"), p_all,
       width = 9.17, height = 3.35)

################################################################################
# 2. Fit both cohorts
#
#    We keep `sigs_fixed = TRUE`
################################################################################
controls <- SignaturePPF_control(maxiter = 4000, tol = 1e-6)

fit_ppf <- function(data, out_file) {
  if (file.exists(out_file)) {
    message("using existing fit: ", basename(out_file))
    return(readRDS(out_file))
  }
  # NOT pruned. The two cohorts are compared reference by reference, and they do
  # NOT switch off the same ones - that difference is the result. Pruning each fit
  # to its own live set would leave the two with different signatures and nothing
  # to line up.
  fit <- SignaturePPF(data,
                      prune_solution = FALSE,
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

p_mu <- ggplot(data = mu_table) +
  geom_abline(color = "gray")+
  geom_point(aes(x = mu_Breast80, y = mu_ICGC, color = signature))+
  scale_color_manual(name = "Signature", values = ALLUVIAL_SIG_COLS)+
  theme_bw() +
  xlab("Relevance weights (Breast80)") +
  ylab("Relevance weights (BreastICGC)")

## The same panel as plot_Mu(), but for both cohorts at once. Built here from
## mu_table rather than by pasting two plot_Mu() calls together with `+`: two
## plots carry two independent scales, so a point of a given size would mean a
## different mu in each half and the columns would not be comparable. One
## ggplot with a facet_grid puts both cohorts on ONE size scale and ONE fill
## scale, and gets the shared legend for free.
mu_long <- data.frame(
  signature = rep(mu_table$signature, 2),
  cohort = factor(rep(c("BreastICGC", "Breast80"), each = nrow(mu_table)),
                  levels = c("BreastICGC", "Breast80")),
  mu = c(mu_table$mu_ICGC, mu_table$mu_Breast80),
  m = c(mu_table$assigned_ICGC, mu_table$assigned_Breast80))
mu_long$signature <- factor(mu_long$signature, levels = mu_table$signature)
mu_long$compressed <- ifelse(mu_long$m == 0 | mu_long$mu < 0.05,
                             "compressed", "not compressed")

p_mu_grid <- ggplot(mu_long, aes(x = cohort, y = 1, size = mu, fill = m,
                                 shape = compressed)) +
  geom_point(colour = "black", stroke = 0.7) +
  facet_wrap(~ signature, ncol = 1, strip.position = "left") +
  scale_x_discrete(position = "top") +
  scale_size(name = expression(mu[k]), range = c(3, 12)) +
  scale_fill_gradientn(
    name = "N. mutations",
    colours = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
  scale_shape_manual(name = "Compressed",
                     values = c("compressed" = 4, "not compressed" = 21)) +
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

ggsave(file.path(FIG_DIR, "Replication_mu_both_cohorts.pdf"), p_mu_grid,
       width = 3.38, height = 7.10)

################################################################################
# 4. Figures
################################################################################
p_burden <- plot_burden_along_genome(list(Breast80 = data80, BreastICGC = dataICGC))
ggsave(file.path(FIG_DIR, "Replication_burden_along_genome.pdf"), p_burden,
       width = 11, height = 3.2)

# --- per-cohort coefficient heatmaps, each beside its relevance-weight column.
#     cap = 1 bounds the colour scale only; the printed number is the estimate.
panel <- function(fit, data) {
  plot_Betas(fit) +
    plot_vector_facets_x(df_assign(fit, data), levs = colnames(fit$Betas)) +
    plot_layout(widths = c(4, 1))
}
p_80 <- panel(fitPPF_80, data80)
p_icgc <- panel(fitPPF_ICGC, dataICGC)

ggsave(file.path(FIG_DIR, "Replication_betas_Breast80.pdf"), p_80, width = 9, height = 6)
ggsave(file.path(FIG_DIR, "Replication_betas_BreastICGC.pdf"), p_icgc, width = 9, height = 6)

# --- the replication figure
p_rep_sig <- plot_beta_replication(fitPPF_80, fitPPF_ICGC,
                               mu_tol = 0.01,
                               colour_by = "signature",
                               label_x = "Breast80", label_y = "BreastICGC",
                               facet_by = "covariate") +
  theme(axis.text = element_text(size= 8)) + theme(aspect.ratio = 1)

ggsave(file.path(FIG_DIR, "Replication_betas_scatter.pdf"), p_rep_sig,
       width = 9.58, height = 6.35)

################################################################################
# 5. Agreement summary
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
  sign_agreement = mean(sign(x) == sign(y)))
write.csv(agreement, file.path(DIR_REPLICATION, "beta_agreement.csv"),
          row.names = FALSE)
print(agreement)

message("done: outputs in ", DIR_REPLICATION, " and ", FIG_DIR)
