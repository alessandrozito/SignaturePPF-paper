################################################################################
# Produces: Figure S9
#
# Stability of covariate effects under a growing covariate set
#
# One train/test split of the GENOME (whole megabases held out, stratified by
# chromosome). On the training bins we fit a sequence of nested models:
#
#     model 0 : no covariate effect  (beta fixed at 0, so the intensity is
#               carried by copy number alone - the position-independent null)
#     model m : PPF with the m covariates chosen by forward selection,
#               m = 1, ..., L   (L = number of covariates)
#
# Forward selection adds, at each step, the covariate most correlated with the
# per-bin residual of the current model on the training bins.
#
# Three outputs:
#   (1) each covariate's coefficient along the sequence -> are the betas stable?
#   (2) alluvial plots of mutation attribution, in- and out-of-sample -> does
#       adding covariates move mutations between signatures?
#   (3) per-patient RMSE of the regional mutation rate against the number of
#       covariates, in- and out-of-sample -> do the covariates predict, or fit?
#
# Usage:  Rscript R/Application_stability_of_covariates.R
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggplot2)
  library(ggalluvial)
  library(patchwork)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()
check_inputs(PATH_ICGC10KB)

## ------------------------------------------------------------------ settings
HOLDOUT_FRAC <- 0.20      # share of 1 Mb regions held out, within each chromosome
REGION_WIDTH <- 1e6       # resolution the mutation rate is scored at
MAXITER <- 4000
TOL <- 1e-6

################################################################################
# 1. Data
################################################################################
dataICGC <- load_cohort(PATH_ICGC10KB)
invisible(SignaturePPF_validate(dataICGC))

covariate_names <- colnames(dataICGC$SignalTrack)
L <- length(covariate_names)
n_bins <- nrow(dataICGC$SignalTrack)

CosmicSigs <- COSMIC_v3.4_SBS96_GRCh37[, SIGS_TO_USE]

bin_of_mut <- bin_of_mutation(dataICGC)
obs <- count_by_bin(bin_of_mut, dataICGC$gr_Mutations$sample, n_bins)

################################################################################
# 2. Random training/test split of the genome
################################################################################
hg19 <- BSgenome.Hsapiens.UCSC.hg19
regions <- tileGenome(seqlengths(hg19)[paste0("chr", c(1:22, "X"))],
                      tilewidth = REGION_WIDTH, cut.last.tile.in.chrom = TRUE)
region_of_bin <- findOverlaps(dataICGC$gr_SignalTrack, regions, select = "first")
chrom_of_bin <- as.character(seqnames(dataICGC$gr_SignalTrack))

set.seed(SEED)
test_regions <- unlist(lapply(split(region_of_bin, chrom_of_bin), function(r) {
  r <- unique(r[!is.na(r)])
  if (length(r) > 1) sample(r, floor(HOLDOUT_FRAC * length(r))) else integer(0)
}), use.names = FALSE)

test_bins <- sort(which(region_of_bin %in% test_regions))
train_bins <- setdiff(seq_len(n_bins), test_bins)

message(sprintf("Split: %s training bins, %s held out (%.1f%%).",
                format(length(train_bins), big.mark = ","),
                format(length(test_bins), big.mark = ","),
                100 * length(test_bins) / n_bins))

saveRDS(list(train_bins = train_bins, test_bins = test_bins, seed = SEED),
        file.path(DIR_STABILITY, "train_test_bins.rds.gzip"), compress = "gzip")

gr_train <- dataICGC$gr_Mutations[bin_of_mut %in% train_bins]
gr_test <- dataICGC$gr_Mutations[bin_of_mut %in% test_bins]

################################################################################
# 3. Fit the sequence
#################################################################################
fit_map <- function(covs, out_file, betas_zero = FALSE) {
  if (file.exists(out_file)) {
    message("using existing fit: ", basename(out_file))
    return(readRDS(out_file))
  }
  train <- subset_bins(dataICGC, train_bins, bin_of_mut, covariates = covs)
  # NOT pruned: the twelve models are compared coefficient by coefficient along
  # the selection path, so every one of them has to report the same signatures.
  fit <- SignaturePPF(
    prune_solution = FALSE,
    train,
    sigs = CosmicSigs,
    sigs_fixed = TRUE,
    method = "map",
    controls = SignaturePPF_control(maxiter = MAXITER, tol = TOL,
                                    update_Betas = !betas_zero),
    init = if (betas_zero) {
      SignaturePPF_init(Betas_start = matrix(0, length(covs), ncol(CosmicSigs)))
    } else {
      SignaturePPF_init()
    },
    seed = SEED,
    verbose = TRUE)
  saveRDS(fit, out_file, compress = "gzip")
  fit
}

fits <- vector("list", L + 1L)
fits[[1]] <- fit_map(covariate_names[1],
                     file.path(DIR_STABILITY, "model_00.rds.gzip"),
                     betas_zero = TRUE)

selected <- character(0)
remaining <- covariate_names
sel_path <- data.frame()

# Residual of the current model on the training bins, summed over patients. The
# null model contributes no positional structure, so step 1 selects against the
# observed counts themselves.
obs_train <- rowSums(obs[train_bins, , drop = FALSE])
residual <- obs_train

for (m in seq_len(L)) {
  cors <- cor(residual,
              dataICGC$SignalTrack[train_bins, remaining, drop = FALSE])[1, ]
  best <- names(which.max(abs(cors)))
  sel_path <- rbind(sel_path,
                    data.frame(model = m, covariate = best,
                               abs_cor = unname(abs(cors[best]))))
  selected <- c(selected, best)
  remaining <- setdiff(remaining, best)

  message(sprintf("\nModel %d/%d: adding '%s' (|cor| = %.3f)",
                  m, L, best, abs(cors[best])))

  fits[[m + 1L]] <- fit_map(selected,
                            file.path(DIR_STABILITY,
                                      sprintf("model_%02d.rds.gzip", m)))

  lambda <- predict_lambda_bins(fits[[m + 1L]], dataICGC, train_bins)
  residual <- obs_train - rowSums(lambda)
  rm(lambda)

  if (!length(remaining)) break
}

saveRDS(sel_path, file.path(DIR_STABILITY, "selection_path.rds.gzip"),
        compress = "gzip")
write.csv(sel_path, file.path(DIR_STABILITY, "selection_path.csv"),
          row.names = FALSE)
print(sel_path)

# Labels shared by all three figures: one per model, in order.
model_labels <- c("No covariates", paste0("+", sel_path$covariate))

################################################################################
# OUTPUT 1: each covariate's coefficient along the sequence
################################################################################
p_betas <- plot_beta_path(fits[-1], covariate_order = sel_path$covariate, mu_cutoff = 0)
ggsave(file.path(FIG_DIR, "Stability_betas_sequence.pdf"), p_betas,
       width = 12.5, height = 4.66)

################################################################################
# OUTPUT 2: attribution flow, in- and out-of-sample
################################################################################
assign_along_path <- function(gr) {
  vapply(fits, function(f) assign_mutations(f, gr), character(length(gr)))
}
A_train <- assign_along_path(gr_train)
A_test <- assign_along_path(gr_test)

river_train <- plot_assignment_alluvial(A_train, model_labels, SIGS_TO_USE) +
  facet_wrap(. ~ "In-sample mutation assignment")
river_test <- plot_assignment_alluvial(A_test, model_labels, SIGS_TO_USE) +
  facet_wrap(. ~ "Held-out mutation assignment")

p_river <- river_train + river_test + plot_layout(guides = "collect")
ggsave(file.path(FIG_DIR, "Stability_riverplots.pdf"), p_river,
       width = 11.33, height = 4.66)

# Total fraction of mutations re-assigned at the end of the track
prop.table(table(A_train[, 1] != A_train[, 12]))
# ---> 17.8% of mutations get reassigned after the introduction of all covariates

################################################################################
# OUTPUT 3: RMSE of the regional mutation rate
#################################################################################
grp_in <- region_of_bin[train_bins];  keep_in <- !is.na(grp_in)
grp_out <- region_of_bin[test_bins];  keep_out <- !is.na(grp_out)
obs_mb_in <- rowSums(rowsum(obs[train_bins[keep_in], , drop = FALSE], grp_in[keep_in]))
obs_mb_out <- rowSums(rowsum(obs[test_bins[keep_out], , drop = FALSE], grp_out[keep_out]))

scored <- lapply(seq_along(fits), function(m) {
  message("scoring model ", m - 1L, "/", L)
  lam_in <- predict_lambda_bins(fits[[m]], dataICGC, train_bins)
  rmse_in <- patient_rmse(lam_in, obs, train_bins, region_of_bin)
  pred_mb_in <- rowSums(rowsum(lam_in[keep_in, , drop = FALSE], grp_in[keep_in]))
  rm(lam_in)
  lam_out <- predict_lambda_bins(fits[[m]], dataICGC, test_bins)
  rmse_out <- patient_rmse(lam_out, obs, test_bins, region_of_bin)
  pred_mb_out <- rowSums(rowsum(lam_out[keep_out, , drop = FALSE], grp_out[keep_out]))
  rm(lam_out)
  list(patient = data.frame(model = m - 1L, patient = names(rmse_in),
                            in_sample = rmse_in, out_sample = rmse_out,
                            row.names = NULL),
       mb = rbind(
         data.frame(model = m - 1L, split = "In-sample",
                    region = as.integer(names(obs_mb_in)),
                    observed = obs_mb_in, predicted = pred_mb_in, row.names = NULL),
         data.frame(model = m - 1L, split = "Held-out",
                    region = as.integer(names(obs_mb_out)),
                    observed = obs_mb_out, predicted = pred_mb_out, row.names = NULL)))
})
rmse <- do.call(rbind, lapply(scored, `[[`, "patient"))
mb <- do.call(rbind, lapply(scored, `[[`, "mb"))
mb$split <- factor(mb$split, levels = c("In-sample", "Held-out"))
mb$residual <- mb$observed - mb$predicted
mb$label <- factor(sprintf("M%d %s", mb$model, model_labels[mb$model + 1L]),
                   levels = sprintf("M%d %s", seq_along(model_labels) - 1L, model_labels))

write.csv(rmse, file.path(DIR_STABILITY, "patient_rmse.csv"), row.names = FALSE)
saveRDS(mb, file.path(DIR_STABILITY, "megabase_predictions.rds.gzip"), compress = "gzip")

cohort_summary <- do.call(rbind, lapply(
  split(mb, list(mb$model, mb$split), drop = TRUE), function(d) {
    data.frame(model = d$model[1], split = d$split[1],
               covariate = model_labels[d$model[1] + 1L],
               n_regions = nrow(d),
               cohort_rmse = sqrt(mean(d$residual^2)),
               cohort_cor = cor(d$predicted, d$observed),
               row.names = NULL)
  }))
cohort_summary <- cohort_summary[order(cohort_summary$split, cohort_summary$model), ]
write.csv(cohort_summary, file.path(DIR_STABILITY, "cohort_regional_summary.csv"),
          row.names = FALSE)
print(cohort_summary)

# Make the plot now
SPLIT_COLS <- c("In-sample" = "#A6CEE3", "Held-out" = "#EB6834")

p_rmse_mb <- ggplot(cohort_summary,
                    aes(factor(model), cohort_rmse,
                        colour = split, group = split)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.8) +
  scale_colour_manual(name = NULL, values = SPLIT_COLS) +
  labs(x = NULL, y = "RMSE per Mb") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank())

p_box <- ggplot(mb, aes(x = factor(model), y = residual, fill = split)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = 2) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = SPLIT_COLS) +
  scale_x_discrete(labels = model_labels) +
  # The colours are named by the legend of the panel directly above, and the two
  # share an x axis, so a second legend saying the same thing is noise.
  guides(fill = "none") +
  labs(x = NULL, y = "Observed - predicted\nmutations per Mb (cohort)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_mb <- p_rmse_mb / p_box +
  plot_layout(heights = c(1, 2.2), guides = "collect") &
  theme(legend.position = "top")
ggsave(file.path(FIG_DIR, "Stability_megabase_residuals.pdf"), p_mb,
       width = 4.66, height = 4.66)

p_rmse <- plot_rmse_path(rmse, model_labels)
ggsave(file.path(FIG_DIR, "Stability_rmse.pdf"), p_rmse, width = 8, height = 5)

# Median across patients, which is what the text quotes.
rmse_summary <- do.call(rbind, lapply(split(rmse, rmse$model), function(d) {
  data.frame(model = d$model[1],
             covariate = model_labels[d$model[1] + 1L],
             in_sample = median(d$in_sample),
             out_sample = median(d$out_sample))
}))
write.csv(rmse_summary, file.path(DIR_STABILITY, "rmse_summary.csv"),
          row.names = FALSE)
print(rmse_summary)


################################################################################
# 4. Save all the objects
################################################################################
saveRDS(list(sel_path = sel_path, rmse = rmse, rmse_summary = rmse_summary,
             mb = mb, cohort_summary = cohort_summary,
             A_train = A_train, A_test = A_test, model_labels = model_labels),
        file.path(DIR_STABILITY, "stability_outputs.rds.gzip"),
        compress = "gzip")

message("done: outputs in ", DIR_STABILITY, " and ", FIG_DIR)
