################################################################################
# Produces: Figure 1
#
# Mutations along the genome
#
#   panel a   figures/Figure1_a_total_mutations_Mb.pdf
#             substitutions per megabase, 113 ICGC Breast-AdenoCa genomes
#   panel b   figures/Figure1_b_2kb_mutations_H3K9me3.pdf
#             chrX:88-89 Mb at 2 kb, over the H3K9me3 track
#
#
# Usage:  Rscript R/Figure1_mutation_landscape.R
#
# chrX:88-89 Mb is the megabase with the strongest correlation between mutation
# count and H3K9me3 among those with at least 200 mutations. That search is
# scan_Mb_correlations() at the bottom; it is slow and does not run on sourcing.
################################################################################

suppressPackageStartupMessages({
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
  library(patchwork)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")

## Checked against its coordinates below, so a change in the tiling fails loudly.
REGION_HIGHLIGHT <- 2986L
REGION_EXPECTED  <- "chrX:88000001-89000000"

BAND_COLS <- c("#4682B4", "antiquewhite2")   # alternating chromosome bands
TRACK_COL <- "#4682B4"

## data/data_for_figure1/ holds just the three files this script needs; fall back
## to the full application set when it is absent.
fig1_input <- function(path) {
  local <- file.path(DIR_DATA_FIGURE1, basename(path))
  if (file.exists(local)) local else path
}
PATH_SNV   <- fig1_input(PATH_ICGC_SNV)
PATH_BLACK <- fig1_input(PATH_BLACKLIST)
PATH_H3K9  <- fig1_input(path_mark("H3K9me3", "tissue"))

check_inputs(c(PATH_SNV, PATH_BLACK, PATH_H3K9))

################################################################################
# 1. Mutations, with the ENCODE blacklist removed
################################################################################
message("== 1. inputs ==")

genome        <- BSgenome.Hsapiens.UCSC.hg19
chrom_lengths <- seqlengths(genome)[1:23]          # chr1-22 and chrX

gr_tumor  <- readRDS(PATH_SNV)
blacklist <- rtracklayer::import(PATH_BLACK)

n_before <- length(gr_tumor)
gr_tumor <- IRanges::subsetByOverlaps(gr_tumor, blacklist, invert = TRUE)
message(sprintf("mutations: %s, %s dropped by the blacklist, %s samples",
                format(n_before, big.mark = ","),
                format(n_before - length(gr_tumor), big.mark = ","),
                length(unique(mcols(gr_tumor)$sample))))

## Same 23 chromosomes, so a 1 Mb tile is exactly 500 of the 2 kb tiles.
hg19_Mb  <- tileGenome(chrom_lengths, tilewidth = 1e6, cut.last.tile.in.chrom = TRUE)
hg19_2kb <- tileGenome(chrom_lengths, tilewidth = 2000, cut.last.tile.in.chrom = TRUE)

region_filter <- hg19_Mb[REGION_HIGHLIGHT]
region_label  <- sprintf("%s:%d-%d", as.character(seqnames(region_filter)),
                         start(region_filter), end(region_filter))
if (!identical(region_label, REGION_EXPECTED)) {
  stop("megabase ", REGION_HIGHLIGHT, " is ", region_label,
       ", not ", REGION_EXPECTED, " - the tiling has changed.", call. = FALSE)
}
message("zoom region: ", region_label)

################################################################################
# 2. Figure 1a - cohort burden per megabase
################################################################################
message("== 2. panel a ==")

df_tumor_all <- as.data.frame(gr_tumor) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor, hg19_Mb))) %>%
  group_by(region) %>%
  summarise(n = n(), .groups = "drop")

## One band per chromosome, alternating, spanning the megabases it contains.
chrom_bands <- as.data.frame(hg19_Mb) %>%
  mutate(region = seq_along(hg19_Mb)) %>%
  group_by(seqnames) %>%
  summarise(xmin = min(region), xmax = max(region) + 1, .groups = "drop") %>%
  arrange(seqnames) %>%
  mutate(fill_color = factor(seq_along(seqnames) %% 2))

panelA <- ggplot() +
  theme_bw(base_size = 11) +
  geom_rect(data = chrom_bands,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = fill_color),
            alpha = 0.2) +
  scale_fill_manual(values = BAND_COLS) +
  geom_point(data = df_tumor_all, aes(x = region, y = n),
             size = 1, shape = 21, color = "gray25") +
  ## The megabase panel b expands, marked in red.
  geom_point(data = df_tumor_all %>% filter(region == REGION_HIGHLIGHT),
             aes(x = region + 0.5, y = n), size = 2, color = "red") +
  xlab("Genomic region (1 Mb)") +
  ylab("Number of mutations") +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0.01, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02)) +
  facet_wrap(. ~ "Total mutations in Breast-AdenoCA (Mb scale)")

## Figure 1, panel a -> img/Figure1_TopoMutSig_breast.pdf (top)
ggsave(file.path(FIG_DIR, "Figure1_a_total_mutations_Mb.pdf"), panelA,
       width = 9.16, height = 2.27)

################################################################################
# 3. Figure 1b - that megabase at 2 kb, against H3K9me3
################################################################################
message("== 3. panel b ==")

signalH3K9me3 <- rtracklayer::import(PATH_H3K9,
                                     which = region_filter)

gr_tumor_filter <- IRanges::subsetByOverlaps(gr_tumor, region_filter)
regions_all     <- subjectHits(findOverlaps(region_filter, hg19_2kb))

## Complete over (sample, bin) before summing, so an unmutated bin is a zero
## rather than an absent row - the open circles in the panel.
df_aggr <- as.data.frame(gr_tumor_filter) %>%
  mutate(region = subjectHits(findOverlaps(gr_tumor_filter, hg19_2kb))) %>%
  group_by(sample, region) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(sample, region = regions_all, fill = list(n = 0)) %>%
  group_by(region) %>%
  summarise(all = sum(n), .groups = "drop") %>%
  mutate(region2 = start(region_filter) + 2000 * seq_along(regions_all) - 2000,
         H3K9me3 = signalH3K9me3$score)

stopifnot(nrow(df_aggr) == length(regions_all),
          length(signalH3K9me3) == length(regions_all))

panelB_top <- df_aggr %>%
  mutate(is_zero = (all == 0)) %>%
  ggplot(aes(x = region2 + 0.5, y = all)) +
  geom_segment(aes(x = region2 + 0.5, xend = region2 + 0.5, y = 0, yend = all),
               linewidth = 0.2, color = TRACK_COL, alpha = 0.7) +
  geom_point(aes(color = is_zero, shape = is_zero), size = 1.2) +
  scale_color_manual(values = c(TRACK_COL, "gray25")) +
  scale_shape_manual(values = c(1, NA)) +     # empty bins draw no point
  theme_bw(base_size = 11) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  ylab("Number of\nmutations") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02))

panelB_bottom <- ggplot(df_aggr) +
  geom_area(aes(x = region2, y = H3K9me3), alpha = 0.9, fill = TRACK_COL) +
  theme_bw(base_size = 11) +
  xlab(sprintf("Genomic region in %s (2 Kb)", as.character(seqnames(region_filter)))) +
  ylab("Average\nfold change") +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.03, 0.02))

panelB <- (panelB_top / panelB_bottom) + plot_layout(heights = c(1, 1))

## Figure 1, panel b -> img/Figure1_TopoMutSig_breast.pdf (bottom)
ggsave(file.path(FIG_DIR, "Figure1_b_2kb_mutations_H3K9me3.pdf"), panelB,
       width = 8.43, height = 2.31)

message("correlation in this megabase: ",
        sprintf("%.3f", cor(df_aggr$all, df_aggr$H3K9me3)))

################################################################################
# 4. How the zoom region was chosen - not run
################################################################################
## Correlation between mutation count and H3K9me3 within each megabase, at 2 kb.
## One bigWig read per megabase, so about an hour. REGION_HIGHLIGHT is its answer.
##
##   corrs <- scan_Mb_correlations()
##   corrs %>% left_join(df_tumor_all, by = "region") %>%
##     filter(n >= 200) %>% arrange(desc(cor_H3K9me3)) %>% head(20)
scan_Mb_correlations <- function(min_mutations = 1L) {
  bw <- PATH_H3K9
  out <- vector("list", length(hg19_Mb))

  for (i in seq_along(hg19_Mb)) {
    reg <- hg19_Mb[i]
    gr_i <- IRanges::subsetByOverlaps(gr_tumor, reg)
    if (length(gr_i) < min_mutations) next

    bins <- subjectHits(findOverlaps(reg, hg19_2kb))
    agg <- as.data.frame(gr_i) %>%
      mutate(region = subjectHits(findOverlaps(gr_i, hg19_2kb))) %>%
      group_by(sample, region) %>%
      summarise(n = n(), .groups = "drop") %>%
      tidyr::complete(sample, region = bins, fill = list(n = 0)) %>%
      group_by(region) %>%
      summarise(all = sum(n), .groups = "drop")

    signal <- rtracklayer::import(bw, which = reg)
    if (length(signal) != nrow(agg)) next     # ragged last tile of a chromosome

    out[[i]] <- data.frame(region = i,
                           cor_H3K9me3 = cor(agg$all, signal$score))
    if (i %% 100 == 0) message("  megabase ", i, " / ", length(hg19_Mb))
  }
  dplyr::bind_rows(out)
}

cat("\nDONE\n")
