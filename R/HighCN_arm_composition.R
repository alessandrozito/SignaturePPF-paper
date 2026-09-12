################################################################################
# Produces: no figure. Writes output/GoodnessOfFit/arm_composition_cells.csv
#
# Arm composition of the volcano: how much of each side sits on an implausible
# consensus copy-number segment.
#
# The negative arm is expected to be dominated by segments the PCAWG consensus
# calls at tens to hundreds of copies. build_CopyTrack() floors copy number at
# 0.1 but does not cap it, so such a segment carries a huge exposure and the
# model predicts mutations that are not there.
#
# Usage:  Rscript R/HighCN_arm_composition.R
################################################################################

suppressPackageStartupMessages({
  library(SignaturePPF)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(tidyverse)
})

## Run from the repository root, or from R/.
source(if (file.exists("config.R")) "config.R" else "../config.R")
load_functions()

WIDTH_MAIN <- 1e6
ALPHA      <- 0.05
CN_CUTS    <- c(50, 100)   # reported side by side; 50 is the headline

data <- readRDS(PATH_ICGC2KB)
fit  <- relabel_by_mu(prune_signatures(
  readRDS(file.path(DIR_DENOVO, "MCMCSolution.rds.gzip"))))

g  <- gof_setup(data, fit, widths = c(1e4, 1e5, WIDTH_MAIN))
rp <- gof_region_patient(g, WIDTH_MAIN)
message("cells: ", nrow(rp))
saveRDS(rp, file.path(DIR_GOF, "highcn_rp_1Mb.rds"), compress = "gzip")

################################################################################
# The high-copy segments, as (sample, 1 Mb window) cells
################################################################################
cn <- readr::read_tsv(file.path(DATA_DIR, "20170119_final_consensus_copynumber_donor"),
                      col_types = readr::cols_only(
                        sampleID = "c", chr = "c", start = "d", end = "d",
                        total_cn = "d"))
cn <- cn[cn$sampleID %in% unique(rp$sample) & !is.na(cn$total_cn), ]

## A window is flagged when any high-copy segment of that donor overlaps it;
## cov_frac records how much of the window the segment actually covers, so a
## one-base clip is distinguishable from a window that is entirely amplified.
flag_cells <- function(cut, cells) {
  hi <- cn[cn$total_cn >= cut, ]
  gr_hi <- GRanges(paste0("chr", hi$chr),
                   IRanges(hi$start, hi$end), sample = hi$sampleID)
  gr_cell <- GRanges(cells$chrom, IRanges(cells$start, width = WIDTH_MAIN))
  ov <- findOverlaps(gr_cell, gr_hi)
  ov <- ov[cells$sample[queryHits(ov)] == gr_hi$sample[subjectHits(ov)]]
  stopifnot(length(ov) > 0)
  w <- as.numeric(width(pintersect(gr_cell[queryHits(ov)], gr_hi[subjectHits(ov)])))
  agg <- tapply(w, queryHits(ov), sum)
  idx <- as.integer(names(agg))
  out <- cells[idx, c("sample", "chrom", "start"), drop = FALSE]
  out$cov_frac <- pmin(as.numeric(agg) / WIDTH_MAIN, 1)
  out
}

################################################################################
# Composition of each arm
################################################################################
report <- function(cut) {
  hi <- flag_cells(cut, unique(rp[, c("sample", "chrom", "start")]))
  d <- rp %>%
    left_join(hi, by = c("sample", "chrom", "start")) %>%
    mutate(high_cn = !is.na(cov_frac),
           sig     = .data$padj < ALPHA & is.finite(.data$padj),
           arm     = ifelse(.data$observed < .data$expected, "negative", "positive"))

  cat(sprintf("\n=================  total_cn >= %d  =================\n", cut))
  cat(sprintf("flagged (sample, Mb) cells: %d of %d (%.3f%%)\n",
              sum(d$high_cn), nrow(d), 100 * mean(d$high_cn)))

  tab <- d %>% filter(sig) %>% group_by(arm) %>%
    summarise(n = n(),
              n_high = sum(high_cn),
              pct_high = 100 * mean(high_cn),
              .groups = "drop")
  cat("\nsignificant cells by arm:\n"); print(as.data.frame(tab), row.names = FALSE)

  ## The converse: of the flagged cells, how many actually show up as
  ## significantly depleted. A high rate is what makes the mechanism credible.
  conv <- d %>% filter(high_cn) %>%
    summarise(n = n(),
              sig_neg = sum(sig & arm == "negative"),
              pct_sig_neg = 100 * mean(sig & arm == "negative"))
  cat("\nof the flagged cells:\n"); print(as.data.frame(conv), row.names = FALSE)

  ## Ranked by evidence, so the claim can be stated for the extreme tail even
  ## if it weakens across the whole arm.
  cat("\nnegative arm, share on a high-CN segment by rank of padj:\n")
  neg <- d %>% filter(sig, arm == "negative") %>% arrange(.data$padj)
  for (k in c(10, 25, 50, 100, 250, nrow(neg))) {
    if (k > nrow(neg)) next
    cat(sprintf("  top %5d: %5.1f%%\n", k, 100 * mean(neg$high_cn[seq_len(k)])))
  }
  invisible(d)
}

for (cut in CN_CUTS) d <- report(cut)

write.csv(d %>% filter(sig) %>%
            dplyr::select(sample, chrom, start, observed, expected,
                          log2fc, padj, arm, high_cn, cov_frac),
          file.path(DIR_GOF, "arm_composition_cells.csv"), row.names = FALSE)
cat("\nDONE\n")
