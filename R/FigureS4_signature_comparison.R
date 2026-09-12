################################################################################
# Produces: Figure S4
#
# De novo signatures under PPF, CompressiveNMF and SignatureAnalyzer
#
# The three panels are matched one-to-one against the PPF solution, so a row is
# the same process across methods; a method with fewer signatures leaves the
# unmatched rows blank. Each row is annotated with its closest COSMIC signature
# and the cosine to it, and the rows matching SBS1 are starred.
#
# Usage:  Rscript R/FigureS4_signature_comparison.R
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF); library(GenomicRanges); library(tidyverse)
  library(patchwork)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()
source(file.path(R_DIR, "Simulation_functions.R"))   # match_MutSign

STAR_SBS <- "SBS1"          # starred in magenta, as in the published figure
STAR_COL <- "magenta3"

################################################################################
# 1. The three solutions
################################################################################
data <- readRDS(PATH_ICGC2KB)
fit  <- relabel_by_mu(prune_signatures(
  readRDS(file.path(DIR_DENOVO, "MCMCSolution.rds.gzip"))))

## The same second permutation the main de novo figures apply, so a row here is
## the row of that figure. Keep in step with Reproduce_figures_Application_denovo.R.
NEW_ORDER <- c(1, 4, 2, 5, 3, 7, 6, 8, 9, 10)
fit <- relabel_by_mu(fit, order = NEW_ORDER)
SigsMean <- fit$Signatures

nmf_dir     <- file.path(DIR_DENOVO, "BaselineNMF")
BaselineNMF <- readRDS(file.path(nmf_dir, "out_CompNMF.rds"))
SigAnalyzer <- readRDS(file.path(nmf_dir, "out_SigAnalyzerL1KL.rds"))

mutMatrix <- getTotalMutations(data$gr_Mutations)

################################################################################
# 2. Reconstruction of the aggregated count matrix
################################################################################
rmse <- c(
  SignaturePPF      = sqrt(mean((fit$Signatures %*% fit$Thetas - mutMatrix)^2)),
  CompressiveNMF    = sqrt(mean((BaselineNMF$Signatures %*% BaselineNMF$Theta -
                                   mutMatrix)^2)),
  SignatureAnalyzer = sqrt(mean((SigAnalyzer$Signature.norm %*%
                                   SigAnalyzer$Exposure - mutMatrix)^2)))
K <- c(SignaturePPF      = ncol(fit$Signatures),
       CompressiveNMF    = ncol(BaselineNMF$Signatures),
       SignatureAnalyzer = ncol(SigAnalyzer$Signature.norm))

summary_tbl <- data.frame(model = names(rmse), K = as.integer(K),
                          rmse = round(as.numeric(rmse), 2))
print(summary_tbl, row.names = FALSE)
write.csv(summary_tbl, file.path(DIR_DENOVO, "count_matrix_rmse.csv"),
          row.names = FALSE)

################################################################################
# 3. Match every method to the PPF rows
################################################################################
## match_MutSign pads the smaller solution with zero columns, which is what
## leaves a row blank rather than shifting the ones below it up.
align <- function(R_hat) {
  m <- match_MutSign(R_true = SigsMean, R_hat = R_hat)
  colnames(m$R_hat) <- colnames(SigsMean)
  m$R_hat
}
panels <- list(
  "Poisson Process Factorization" = SigsMean,
  "CompNMF"                       = align(BaselineNMF$Signatures),
  "SignatureAnalyzer"             = align(SigAnalyzer$Signature.norm))

## Rows carrying no signature under a given method keep their slot but lose their
## number, so the three panels stay row-aligned. Empty rows get a unique run of
## spaces: facet levels must be distinct, but these render as blank strips.
row_labels <- function(M) {
  filled <- colSums(M) > 0
  lab <- strrep(" ", seq_len(ncol(M)))
  lab[filled] <- as.character(seq_len(sum(filled)))
  lab
}

PALETTE <- c("#40BDEE", "#020202", "#E52925", "#CCC9CA", "#A3CF62", "#ECC5C5")

################################################################################
# 4. One panel
################################################################################
build_panel <- function(M, title) {
  filled <- colSums(M) > 0
  lab <- row_labels(M)
  Mp <- M
  colnames(Mp) <- lab

  ann <- match_to_cosmic(M[, filled, drop = FALSE])
  ann$Signature <- factor(lab[filled], levels = lab)
  ann$Mutation  <- factor("T>G", levels = c("C>A", "C>G", "C>T",
                                            "T>A", "T>C", "T>G"))
  ann$starred <- ann$best_match == STAR_SBS
  ann$text <- ifelse(ann$starred,
                     sprintf("* %s (%.2f)", ann$best_match, ann$cosine),
                     sprintf("%s (%.2f)",   ann$best_match, ann$cosine))

  p <- plot_Signatures(Mp, ci = FALSE) +
    ## coord_cartesian(clip = "off") keeps the longer COSMIC names (SBS40a,
    ## SBS10b) from being cut at the panel edge.
    coord_cartesian(clip = "off") +
    geom_text(data = ann, inherit.aes = FALSE,
              aes(x = Inf, y = Inf, label = text, colour = starred),
              fontface = ifelse(ann$starred, 2, 1),
              hjust = 1, vjust = 1.3, size = 2.9) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = STAR_COL),
                        guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12,
                                    colour = "grey25"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, colour = "grey45"),
          legend.position = "none",
          plot.margin = margin(4, 18, 10, 4))

  ## The top strips carry the mutation class, coloured as the bars are.
  if (requireNamespace("ggtext", quietly = TRUE)) {
    mut <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    labs_md <- stats::setNames(
      sprintf("<span style='color:%s'>**%s**</span>", PALETTE,
              gsub(">", "&gt;", mut)), mut)
    p <- p +
      facet_grid(Signature ~ Mutation, scales = "free", switch = "y",
                 labeller = labeller(Mutation = labs_md)) +
      theme(strip.text.x = ggtext::element_markdown(size = 9))
  }
  p
}

p <- build_panel(panels[[1]], names(panels)[1]) |
     build_panel(panels[[2]], names(panels)[2]) |
     build_panel(panels[[3]], names(panels)[3])

ggsave(file.path(FIG_DIR, "Breast_suppl_Signatures_comparison.pdf"), plot = p,
       width = 13.42, height = 7.3)
message("written: Breast_suppl_Signatures_comparison.pdf")
