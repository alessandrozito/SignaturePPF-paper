################################################################################
# Produces: Figure 2, and Figure S7
#
# Goodness of fit of the PPF model
#
#   top     Pearson residuals at 1 Mb, along the genome
#   bottom  observed vs predicted with and without covariates | (patient, Mb)
#           volcano | KS plot for the five most deviating patients
#
# Usage:  Rscript R/Figure2_goodness_of_fit.R
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(tidyverse)
  library(patchwork)
  library(ggnewscale)
})

source(file.path(Sys.getenv("SIGNATUREPPF_PAPER",
                            unset = path.expand("~/SignaturePPF-paper")),
                 "config.R"))
load_functions()

WIDTH_MAIN <- 1e6
ALPHA <- 0.05
DRIVER_CUT <- 0.5          # share of a region's excess above which one patient owns it
N_LABEL_VOLCANO <- 8       # cells named in the volcano
N_KS <- 5                  # worst patients drawn in the uniform panel

cache <- function(file, expr, require_cols = NULL) {
  path <- file.path(DIR_GOF, file)
  if (file.exists(path)) {
    out <- readRDS(path)
    ok <- is.null(require_cols) ||
      (is.data.frame(out) && all(require_cols %in% names(out)))
    if (ok) { message("using cached ", file); return(out) }
    message("cache ", file, " is stale (missing: ",
            paste(setdiff(require_cols, names(out)), collapse = ", "),
            ") - recomputing")
  }
  out <- force(expr); saveRDS(out, path, compress = "gzip"); out
}

################################################################################
# 1. Data and fits
################################################################################
data <- readRDS(PATH_ICGC2KB)
fit <- relabel_by_mu(prune_signatures(
  readRDS(file.path(DIR_DENOVO, "MCMCSolution.rds.gzip"))))
outNoCovs <- relabel_by_mu(readRDS(file.path(DIR_DENOVO,
                                             "MAPSolution_noCovariates.rds.gzip")))

g <- gof_setup(data, fit, widths = c(1e4, 1e5, WIDTH_MAIN))

################################################################################
# 2. Observed against predicted, at the megabase scale
################################################################################
hg19_Mb <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23],
                      tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
over <- findOverlaps(data$gr_SignalTrack, hg19_Mb)
df_chrom <- as.data.frame(hg19_Mb) %>% mutate(region = seq_len(n()))

agg_pred <- function(f) {
  tr <- data$gr_SignalTrack
  tr$LambdaPred <- rowSums(reconstruct_lambda(f, data$SignalTrack, data$CopyTrack))
  as.data.frame(tr) %>%
    mutate(region = subjectHits(over)) %>%
    group_by(region) %>%
    summarize(Lambda = sum(LambdaPred), .groups = "drop")
}

count_windows <- function(df_lambda) {
  as.data.frame(data$gr_Mutations) %>%
    mutate(region = subjectHits(findOverlaps(data$gr_Mutations, hg19_Mb))) %>%
    group_by(region) %>%
    summarise(n = n(), .groups = "drop") %>%
    left_join(df_lambda, by = "region") %>%
    left_join(dplyr::select(df_chrom, region, seqnames), by = "region") %>%
    mutate(chr_num = as.integer(gsub("chr", "", seqnames)),
           point_color_group = as.factor(ifelse(is.na(chr_num %% 2), 1,
                                                chr_num %% 2)))
}

df_ppf <- count_windows(agg_pred(fit))
df_nocov <- count_windows(agg_pred(outNoCovs))

chrom_bands <- as.data.frame(hg19_Mb) %>%
  mutate(region = seq_along(hg19_Mb)) %>%
  group_by(seqnames) %>%
  summarise(xmin = min(region), xmax = max(region) + 1,
            mid = mean(range(region)), .groups = "drop") %>%
  arrange(seqnames) %>%
  mutate(fill_color = factor(seq_along(seqnames) %% 2),
         label = sub("^chr", "", as.character(seqnames)))

## The covariate-free prediction is drawn underneath in grey, so the comparison
## is inside the panel rather than a second row of the figure.
plot_track_both <- function(df, df_nocov) {
  ggplot() +
    geom_rect(data = chrom_bands,
              aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf,
                  fill = fill_color), alpha = 0.1) +
    scale_fill_manual(values = c("#4682B4", "antiquewhite"), guide = "none") +
    geom_point(data = df, size = 0.66,
               aes(x = region, y = n), shape = 1, color = "gray40") +
    #scale_color_manual(values = c("#93BAF1", "#93BAF1"), guide = "none") +
    ggnewscale::new_scale_color() +
    geom_line(data = df, aes(x = region, y = Lambda, group = seqnames,
                             color = "PPF, covariates + copy num."),
              linewidth = 0.8, alpha = 0.5) +
    geom_line(data = df_nocov, aes(x = region, y = Lambda, group = seqnames,
                                   color = "PPF, copy num. only"),
              linewidth = 0.55, alpha = 0.95) +
    scale_color_manual(values = c("PPF, covariates + copy num." = "tomato",
                                  "PPF, copy num. only" = "darkred"), name = NULL) +
    scale_x_continuous(expand = c(0.01, 0.02), breaks = chrom_bands$mid,
                       labels = chrom_bands$label) +
    scale_y_continuous(expand = c(0.03, 0.02)) +
    labs(y = "Number of mutations", x = "Chromosome") +
    theme_bw()+
    theme(legend.position = "top", legend.text = element_text(size = 10))
    #theme(legend.position = c(0.99, 0.97), legend.justification = c(1, 1),
    #      legend.background = element_rect(fill = "white", colour = "grey70"),
    #      legend.key.size = unit(12, "pt"), legend.text = element_text(size = 8.5),
    #      panel.grid.minor = element_blank(),
    #      panel.grid.major.x = element_blank())
}



plot_reconstruction <- function(df, lim = c(0, 610)) {
  ggplot(df) +
    geom_point(aes(x = Lambda, y = n, color = point_color_group),
               alpha = 0.25, size = 0.5, shape = 20) +
    scale_color_manual(values = c("#000D8B", "#93BAF1"), guide = "none") +
    geom_abline(slope = 1, intercept = 0, color = "#CD2626", linewidth = 0.7) +
    xlim(lim) + ylim(lim) +
    labs(y = "Observed mutations", x = "Predicted mutations") +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid.minor = element_blank())
}

# Plot all three panels
(p_Track <- plot_track_both(df_ppf, df_nocov))
(p_reconstr <- plot_reconstruction(df_ppf))
(p_reconstrNoCovs <- plot_reconstruction(df_nocov))

p_Track + p_reconstr + p_reconstrNoCovs
ggsave(file.path(FIG_DIR, "Figure3_recontructed_Mbscale_join.pdf"),
       width = 14.99, height = 3.71)


################################################################################
# 3. Residuals along the genome
################################################################################
## NOT `res`: that name is taken by the patient results further down, and
## everything below reads `reg`.
reg <- gof_regions(g, width = WIDTH_MAIN)

rp <- gof_region_patient(g, WIDTH_MAIN)
drv <- gof_region_driver(rp)

lab_file <- file.path(DIR_GOF, "region_labels_1Mb.csv")
if (file.exists(lab_file)) {
  lt <- utils::read.csv(lab_file, stringsAsFactors = FALSE)
  labels <- lt$label; classes <- lt$class
} else {
  top_idx <- order(-reg$contrib)[seq_len(min(60L, nrow(reg)))]
  classes <- labels <- rep("", nrow(reg))
  classes[top_idx] <- annotate_regions_class(reg[top_idx, ])
  labels[top_idx] <- region_labels(reg[top_idx, ])
  utils::write.csv(data.frame(region = reg$region, label = labels,
                              class = classes), lab_file, row.names = FALSE)
}

tt <- order(-reg$contrib)[seq_len(min(15L, nrow(reg)))]
message(sprintf(
  "\n15 worst megabases: %d in an immunoglobulin locus; usable sequence %.1f-%.1f%% (cohort median %.1f%%)",
  sum(classes[tt] %in% c("IGH", "IGK", "IGL")),
  100 * min(reg$usable_frac[tt]), 100 * max(reg$usable_frac[tt]),
  100 * stats::median(reg$usable_frac)))

## Regions whose excess is one patient's are drawn open rather than filled.
sig <- !is.na(reg$padj) & reg$padj < ALPHA
sh <- drv$driver_share[match(reg$region, drv$region)]
one_pt <- sig & reg$resid > 0 & !is.na(sh) & sh > DRIVER_CUT
message(sprintf(
  "%d of %d Mb windows rejected; %d of the %d positive ones are one patient's",
  sum(sig), nrow(reg), sum(one_pt), sum(sig & reg$resid > 0)))

drv_tbl <- drv[match(reg$region[one_pt], drv$region), ]
drv_tbl$chrom <- reg$chrom[one_pt]; drv_tbl$start <- reg$start[one_pt]
drv_tbl$resid <- reg$resid[one_pt]
utils::write.csv(drv_tbl[order(-drv_tbl$driver_excess), ],
                 file.path(DIR_GOF, "regions_driven_by_one_patient.csv"),
                 row.names = FALSE)

top <- reg[order(-reg$contrib), ][seq_len(min(25L, nrow(reg))), ]
top$class <- classes[top$region]
utils::write.csv(
  top[, c("chrom", "start", "end", "observed", "expected", "resid",
          "contrib_share", "padj", "usable_frac", "class")],
  file.path(DIR_GOF, "top_misfit_regions.csv"), row.names = FALSE)
print(utils::head(top[, c("chrom", "start", "observed", "expected", "resid",
                          "usable_frac", "class")], 10),
      row.names = FALSE, digits = 3)

ggsave(file.path(FIG_DIR, "FigureS_residuals_along_genome.pdf"),
       plot_gof_residual_track(reg, alpha = ALPHA, labels = labels,
                               driver = drv, driver_cut = DRIVER_CUT),
       width = 13, height = 4.6)

################################################################################
# 4. Per-patient tests
################################################################################
pf <- file.path(DIR_GOF, "fig3_patients.rds")
if (file.exists(pf)) {
  res <- readRDS(pf)
} else {
  res <- gof_patients(g, width = WIDTH_MAIN, n0 = NULL, n_rep = 25,
                      seed = SEED, qq = TRUE)
  saveRDS(res, pf, compress = "gzip")
}
pat <- res$patients
utils::write.csv(pat, file.path(DIR_GOF, "patient_deviations.csv"),
                 row.names = FALSE)
message(sprintf(
  "patients rejected: %d of %d on all their mutations | %d of %d after thinning to n0 = %d",
  sum(pat$padj_full < ALPHA, na.rm = TRUE), nrow(pat),
  sum(pat$padj_thin < ALPHA, na.rm = TRUE), nrow(pat),
  max(pat$n_thin, na.rm = TRUE)))

qq_pat <- cache("fig3_qq_patient_full.rds", gof_qq_curves(g, by = "patient"),
                require_cols = "ks_empirical")
qq_all <- cache("fig3_qq_all_full.rds",
                gof_qq_curves(g, by = "all", verbose = FALSE),
                require_cols = "ks_empirical")
message(sprintf("cohort as one process: n = %s, KS distance %.4f (critical %.4f)",
                format(qq_all$n[1], big.mark = ","), qq_all$ks_full[1],
                ks_critical(qq_all$n[1])))

################################################################################
# 5. Volcano
################################################################################
(p_volcano <- plot_gof_volcano(rp, alpha = 0.05, n_label = 10, label_floor = -1.8))

## The five worst patients, ranked by KS distance RELATIVE TO the critical value
## at their own burden.
pat$ks_ratio_crit <- pat$ks_gap / pat$ks_crit
worst <- utils::head(pat$sample[order(-pat$ks_ratio_crit)], N_KS)
message("worst patients: ", paste(worst, collapse = ", "))

## Exponential QQ, one curve per chromosome, rescaled within the chromosome.
qq_chr <- cache("fig3_qq_chrom_full.rds",
                gof_qq_curves(g, by = "chromosome", verbose = FALSE),
                require_cols = "ks_gap")
kc <- unique(qq_chr[, c("group", "n", "ks_gap")])
print(utils::head(kc[order(-kc$ks_gap), ], 5), row.names = FALSE, digits = 3)

(p_exp_chr <- plot_gof_qq_curves(qq_chr, scale = "detrended", n_highlight = 5) +
  theme(aspect.ratio = 1))

## Uniform QQ: the cohort total and the five worst patients, nothing else.
(p_ks <- plot_gof_qq_curves(qq_pat, scale = "detrended", highlight = worst,
                           only_highlight = FALSE, band = TRUE, overlay_colour = "darkorange",
                           overlay = qq_all) +
  theme(aspect.ratio = 1))


## No expected-count filter: the table lists what the volcano draws, and the
## largest departures sit in cells with a small expected count by construction.
vt <- rp[order(rp$p_val), ][seq_len(min(30L, nrow(rp))), ]
vt$class <- classes[vt$region]
utils::write.csv(vt[, c("sample", "chrom", "start", "end", "observed",
                        "expected", "log2fc", "p_val", "padj", "class")],
                 file.path(DIR_GOF, "top_patient_region_cells.csv"),
                 row.names = FALSE)
print(utils::head(vt[, c("sample", "chrom", "start", "observed", "expected",
                         "log2fc", "padj")], 8), row.names = FALSE, digits = 3)

p_volcano + p_ks + p_exp_chr
ggsave(file.path(FIG_DIR, "Figure3_misspecification.pdf"),
       width = 15.43, height = 4.79)


################################################################################
# 6. Supplements
################################################################################
message("\n== residuals by mutation class ==")
reg_class <- cache("fig3_regions_1Mb_by_class.rds", require_cols = "usable_frac",
                   expr = do.call(rbind, lapply(MUT_CLASSES, function(cl)
                     gof_regions(g, WIDTH_MAIN, class = cl))))
drv_class <- cache("fig3_region_driver_by_class.rds", expr = {
  setNames(lapply(MUT_CLASSES, function(cl) {
    message("  ", cl)
    gof_region_driver(gof_region_patient(g, WIDTH_MAIN, class = cl,
                                         compute_p = FALSE))
  }), MUT_CLASSES)
})

one_pt_class <- do.call(rbind, lapply(MUT_CLASSES, function(cl) {
  d <- reg_class[reg_class$class == cl, ]
  s <- !is.na(d$padj) & d$padj < ALPHA & d$resid > 0
  sh <- drv_class[[cl]]$driver_share[match(d$region, drv_class[[cl]]$region)]
  data.frame(class = cl, n_flagged_pos = sum(s),
             n_one_patient = sum(s & !is.na(sh) & sh > DRIVER_CUT),
             median_driver_share = stats::median(sh[s], na.rm = TRUE),
             row.names = NULL)
}))
utils::write.csv(one_pt_class,
                 file.path(DIR_GOF, "one_patient_regions_by_class.csv"),
                 row.names = FALSE)
print(one_pt_class, row.names = FALSE, digits = 3)

message("\n== 2 kb count distribution ==")
counts <- cache("fig3_count_distribution_2kb.rds",
                require_cols = c("ratio", "mut_share"),
                expr = do.call(rbind, lapply(c("all", MUT_CLASSES),
                                             function(cl)
                                               gof_count_distribution(g, cl))))
utils::write.csv(counts, file.path(DIR_GOF, "count_distribution_2kb.csv"),
                 row.names = FALSE)
ca <- counts[counts$class == "all", ]
print(ca[, c("count", "observed", "expected", "ratio", "mutations",
             "mut_share")], row.names = FALSE, digits = 4)
message(sprintf(
  "cells with 0 or 1 mutation hold %.1f%% of the cohort, reproduced to within %.1f%%; the excess grows from %.1f-fold at 2 mutations to %.0f-fold at %s, over %.1f%% of the cohort",
  100 * sum(ca$mut_share[ca$count %in% c("0", "1")]),
  100 * max(abs(1 - ca$ratio[ca$count %in% c("0", "1")])),
  ca$ratio[3], ca$ratio[nrow(ca)], ca$count[nrow(ca)],
  100 * sum(ca$mut_share[-(1:2)])))

p_counts <- plot_gof_counts(counts)
ggsave(file.path(FIG_DIR, "GoF_count_distribution_2kb.pdf"),
       p_counts, width = 6.5, height = 4.5)


p_track_class <- plot_gof_residual_track_class(reg_class, alpha = ALPHA, labels = labels,
                                                     driver = drv_class,
                                                     driver_cut = DRIVER_CUT, ncol = 2)


p_track_class + (p_counts + theme(legend.position = "right")) + plot_layout(widths = c(2,1))
ggsave(file.path(FIG_DIR, "FigureS_residuals_by_mutation_class.pdf"),
       width = 14.60, height = 4.88)

message("\ndone")
