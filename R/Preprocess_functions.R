## Building a binned cohort object from the raw tracks.
##
## Carried over from the predecessor project's loader. One correctness fix, at
## `merge_with_tumor()` - see the note there.

#' Usable sequence in each bin, after assembly gaps and the blacklist
#'
#' A bin is not worth its nominal width: assembly gaps are runs of N, and
#' blacklisted regions are unmappable. Both are removed from the bin's weight,
#' and a bin left with nothing is marked by a weight of zero so the caller can
#' drop it. That weight is the exposure the Poisson process integrates over, so
#' getting it wrong biases every intensity.
add_bin_weights <- function(gr, blacklist_file = PATH_BLACKLIST,
                            gap_file = PATH_GAPS) {
  blacklist <- rtracklayer::import(blacklist_file)
  gap_gr <- rtracklayer::import(gap_file)

  BinWeight <- GenomicRanges::width(gr)
  score <- gr$score

  # Assembly gaps: count the real bases rather than assuming the whole overlap
  # is N, since a gap record can extend past the bin.
  overlapGaps <- GenomicRanges::findOverlaps(gr, gap_gr)
  rangesGaps <- gr[S4Vectors::queryHits(overlapGaps)]
  seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                           rangesGaps)
  countACGT <- rowSums(Biostrings::alphabetFrequency(seqs)[, c("A", "C", "G", "T")])
  BinWeight[S4Vectors::queryHits(overlapGaps)] <- countACGT

  # Blacklist: subtract the intersected width.
  overlap_blackList <- GenomicRanges::findOverlaps(gr, blacklist)
  pairs <- S4Vectors::Pairs(gr[S4Vectors::queryHits(overlap_blackList)],
                            blacklist[S4Vectors::subjectHits(overlap_blackList)])
  grIntersect <- IRanges::pintersect(pairs)
  BinWeight[S4Vectors::queryHits(overlap_blackList)] <-
    pmax(BinWeight[S4Vectors::queryHits(overlap_blackList)] -
           GenomicRanges::width(grIntersect), 0)

  score[BinWeight == 0] <- 0
  gr$score <- score
  gr$bin_weight <- BinWeight
  gr
}


#' Average a scored track onto a set of bins
bin_gr <- function(gr, genome_aggreg, std_chrs) {
  gr <- GenomeInfoDb::keepSeqlevels(gr, std_chrs, pruning.mode = "coarse")
  GenomicRanges::binnedAverage(
    bins = GenomeInfoDb::keepSeqlevels(genome_aggreg, std_chrs,
                                       pruning.mode = "coarse"),
    numvar = GenomicRanges::coverage(gr, weight = "score"),
    varname = "score")
}


#' The genome tiled at a given width, with usable-sequence weights
tiled_genome <- function(tilewidth) {
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(genome)[1:23],
                            tilewidth = tilewidth,
                            cut.last.tile.in.chrom = TRUE)
}


#' Per-patient copy number on the binned genome
#'
#' Halved, so the entry is copies relative to the diploid genome: the intensity
#' carries a factor of one half and this folds it in. Bins below 0.1 copies are
#' floored, since a zero would make the sample contribute no exposure at all
#' there and the copy calls are not that precise.
build_CopyTrack <- function(gr_tumor, gr_copy, tilewidth = 2000) {
  genome_aggreg <- tiled_genome(tilewidth)
  std_chrs <- paste0("chr", c(1:22, "X"))
  gr_CopyTrack <- add_bin_weights(genome_aggreg)
  bin_weight <- gr_CopyTrack$bin_weight

  samples <- unique(gr_tumor$sample)
  gr_copy_tumor <- gr_copy[gr_copy$sample %in% samples]
  gr_copy_tumor$score[is.na(gr_copy_tumor$score)] <- 2

  CopyTrack <- matrix(0, nrow = length(genome_aggreg), ncol = length(samples))
  for (s in seq_along(samples)) {
    gr_copy_sample <- bin_gr(gr_copy_tumor[gr_copy_tumor$sample == samples[s]],
                             genome_aggreg, std_chrs)
    score <- gr_copy_sample$score
    score[score < 0.1] <- 0.1
    CopyTrack[, s] <- score / 2
  }
  colnames(CopyTrack) <- samples
  GenomicRanges::mcols(gr_CopyTrack) <- cbind("bin_weight" = bin_weight, CopyTrack)
  gr_CopyTrack[bin_weight > 0]
}


#' The eleven genomic covariates on the binned genome
#'
#' The seven chromatin marks are the average of a tissue and a cell-line assay;
#' the other four have a single source. Returned unstandardised.
build_SignalTrack <- function(tilewidth = 2000, verbose = TRUE) {
  genome_aggreg <- tiled_genome(tilewidth)
  std_chrs <- paste0("chr", c(1:22, "X"))
  gr_SignalTrack <- add_bin_weights(genome_aggreg)

  say <- function(...) if (verbose) message("  ", ...)

  say("GC")
  gr_SignalTrack$GC <- bin_gr(rtracklayer::import(PATH_GC),
                              genome_aggreg, std_chrs)$score
  say("Methyl")
  gr_SignalTrack$Methyl <- bin_gr(rtracklayer::import(PATH_METHYL),
                                  genome_aggreg, std_chrs)$score

  for (mark in CHROMATIN_MARKS) {
    say(mark)
    tissue <- bin_gr(rtracklayer::import(path_mark(mark, "tissue")),
                     genome_aggreg, std_chrs)$score
    cell <- bin_gr(rtracklayer::import(path_mark(mark, "cell")),
                   genome_aggreg, std_chrs)$score
    GenomicRanges::mcols(gr_SignalTrack)[[mark]] <- (tissue + cell) / 2
  }

  say("RepliTime")
  gr_SignalTrack$RepliTime <- bin_gr(rtracklayer::import(PATH_REPLITIME),
                                     genome_aggreg, std_chrs)$score
  say("NuclOccup")
  gr_SignalTrack$NuclOccup <- bin_gr(rtracklayer::import(PATH_NUCLEOSOME),
                                     genome_aggreg, std_chrs)$score

  # Drop `score`, which add_bin_weights() leaves behind and which is not a
  # covariate; bin_weight stays as the first column.
  keep <- setdiff(names(GenomicRanges::mcols(gr_SignalTrack)), "score")
  GenomicRanges::mcols(gr_SignalTrack) <- GenomicRanges::mcols(gr_SignalTrack)[, keep]
  gr_SignalTrack[gr_SignalTrack$bin_weight > 0]
}


#' Winsorise then standardise a covariate
#'
#' The bigWig tracks have long right tails - a handful of bins carry signal
#' orders of magnitude above the rest - and an unwinsorised covariate lets those
#' few bins set the scale of the coefficient. Capping at the 0.1% quantiles
#' bounds their leverage without discarding them.
standardize_covariates <- function(x, q_min = 0.001, q_max = 0.999) {
  y <- pmin(x, stats::quantile(x, q_max))
  y <- pmax(y, stats::quantile(x, q_min))
  scale(y)
}


#' Attach each mutation's covariate values, dropping mutations off the grid
#'
#' FIXED relative to the predecessor, which pre-filled the covariate matrix with
#' ZEROS and wrote only the rows that matched a bin. The covariates are
#' standardised, so a zero row is not "missing" - it reads as a perfectly
#' average bin, and those mutations were silently fitted as if they sat in one.
#' They are mutations in assembly gaps and blacklisted regions, whose bins were
#' dropped from the SignalTrack precisely because they carry no usable sequence.
#' Here they are dropped from the mutations too, which is the only consistent
#' choice: the model integrates its intensity over the retained bins, so a
#' mutation outside them has no exposure behind it.
merge_with_tumor <- function(gr_tumor, gr_SignalTrack, verbose = TRUE) {
  hit <- GenomicRanges::findOverlaps(gr_tumor, gr_SignalTrack, select = "first")
  off_grid <- is.na(hit)
  if (any(off_grid) && verbose) {
    message(sprintf("  dropping %s of %s mutations that fall outside every ",
                    format(sum(off_grid), big.mark = ","),
                    format(length(gr_tumor), big.mark = ",")),
            "retained bin (assembly gaps and blacklisted regions).")
  }
  gr_tumor <- gr_tumor[!off_grid]
  hit <- hit[!off_grid]

  covariates <- setdiff(names(GenomicRanges::mcols(gr_SignalTrack)), "bin_weight")
  Xmat <- as.matrix(GenomicRanges::mcols(gr_SignalTrack)[hit, covariates, drop = FALSE])
  GenomicRanges::mcols(gr_tumor) <- cbind(GenomicRanges::mcols(gr_tumor), Xmat)
  gr_tumor
}


################################################################################
# The 80-cancer cohort (Davies et al. 2017)
#
# A different front end from the ICGC one: mutations arrive as per-sample CaVEMan
# VCFs rather than an assembled GRanges, and copy number as per-sample ASCAT
# segment tables rather than one consensus file. Everything downstream of
# `build_breast80_dataset()` is the same code the ICGC path uses.
################################################################################

#' The 96 pyrimidine-centred substitution channels, in COSMIC order
mutation_channels_96 <- function() {
  nucleotides <- c("A", "C", "G", "T")
  substitutions <- c("C>A", "C>T", "C>G", "T>A", "T>G", "T>C")
  sort(apply(expand.grid(nucleotides, substitutions, nucleotides), 1,
             function(x) paste0(x[1], "[", x[2], "]", x[3])))
}


#' Assign each SNV its trinucleotide channel
#'
#' The channel is defined on the pyrimidine of the pair, so a mutation with a
#' purine reference is reported on the opposite strand: the context is reverse
#' complemented and the alleles complemented. Without that the same physical
#' event would land in two different channels depending on which strand the
#' caller happened to report.
#'
#' @param gr A `GRanges` of single-base substitutions with `ref` and `alt`.
#' @param genome The reference the contexts are read from.
#' @return `gr` with a `channel` factor over the 96 channels.
call_mutation_channel <- function(
    gr, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19) {
  ctx_gr <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(gr),
    ranges = IRanges::IRanges(start = GenomicRanges::start(gr) - 1, width = 3))
  context <- BSgenome::getSeq(genome, ctx_gr)
  rccontext <- Biostrings::reverseComplement(context)

  ref <- Biostrings::DNAStringSet(gr$ref)
  alt <- Biostrings::DNAStringSet(gr$alt)
  pyrimidine <- as.character(ref) %in% c("C", "T")

  gr$channel <- factor(
    ifelse(pyrimidine,
           paste0(XVector::subseq(context, 1, 1), "[", ref, ">", alt, "]",
                  XVector::subseq(context, 3, 3)),
           paste0(XVector::subseq(rccontext, 1, 1), "[",
                  Biostrings::complement(ref), ">", Biostrings::complement(alt),
                  "]", XVector::subseq(rccontext, 3, 3))),
    levels = mutation_channels_96())
  gr
}


#' Read a directory of CaVEMan VCFs into one GRanges of SNVs
#'
#' Only clean single-base substitutions are kept: an indel or a multi-allelic
#' record has no trinucleotide channel, so it cannot enter the model.
read_caveman_vcfs <- function(dir = PATH_BREAST80_SNV, tumor = "Breast80",
                              verbose = TRUE) {
  files <- list.files(dir, pattern = "\\.caveman\\.vcf$", full.names = TRUE)
  if (!length(files)) stop("no *.caveman.vcf under ", dir, call. = FALSE)
  if (verbose) message("  reading ", length(files), " VCFs")

  bases <- c("A", "C", "G", "T")
  muts <- data.table::rbindlist(lapply(files, function(f) {
    v <- data.table::fread(f, skip = "#CHROM", showProgress = FALSE)
    data.table::setnames(v, "#CHROM", "CHROM")
    # Base `%in%` and plain vector indexing: data.table is used here for fread's
    # speed only, and its non-standard evaluation needs the package attached,
    # which a sourced helper cannot assume.
    keep <- v$REF %in% bases & v$ALT %in% bases
    data.frame(sample = sub("\\.caveman\\.vcf$", "", basename(f)),
               chrom = paste0("chr", v$CHROM[keep]),
               pos = v$POS[keep], ref = v$REF[keep], alt = v$ALT[keep],
               stringsAsFactors = FALSE)
  }))
  muts <- muts[muts$chrom != "chrY", ]

  gr <- GenomicRanges::GRanges(
    seqnames = muts$chrom,
    ranges = IRanges::IRanges(start = muts$pos, width = 1),
    strand = "*",
    tumor = tumor,
    sample = muts$sample,
    ref = muts$ref,
    alt = muts$alt)
  gr <- call_mutation_channel(gr)
  # `ref` and `alt` have done their job. Left in place they would be two extra
  # mcols the model has to be told to ignore.
  gr$ref <- NULL
  gr$alt <- NULL
  if (verbose) {
    message("  ", format(length(gr), big.mark = ","), " SNVs, ",
            length(unique(gr$sample)), " samples")
  }
  gr
}


#' Read a directory of ASCAT segment tables into a copy-number GRanges
#'
#' The files are headerless, and their chromosome column is numeric in some
#' releases and already named in others, so 23/24 are mapped to X/Y explicitly.
read_ascat_segments <- function(dir = PATH_BREAST80_CN, verbose = TRUE) {
  files <- list.files(dir, pattern = "ascat.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("no ascat *.csv under ", dir, call. = FALSE)
  if (verbose) message("  reading ", length(files), " ASCAT tables")

  cols <- c("seg_id", "chr_num", "start", "end", "normal_total", "normal_minor",
            "tumour_total", "tumour_minor")
  df <- do.call(rbind, lapply(files, function(f) {
    d <- utils::read.csv(f, header = FALSE, col.names = cols,
                         colClasses = "character")
    d$sampleID <- sub("[._]ascat.*$", "", basename(f))
    d
  }))

  chr <- as.character(df$chr_num)
  chr[chr == "23"] <- "X"
  chr[chr == "24"] <- "Y"

  GenomicRanges::GRanges(
    seqnames = paste0("chr", chr),
    ranges = IRanges::IRanges(start = as.numeric(df$start),
                              end = as.numeric(df$end)),
    strand = "*",
    sample = df$sampleID,
    score = as.numeric(df$tumour_total))
}


#' Build the binned 80-cancer cohort object
#'
#' @param tilewidth Bin width in bases. 10 kb by default, which is what the
#'   replication analysis needs - it compares this cohort against ICGC on one
#'   grid, so the two have to be binned identically.
build_breast80_dataset <- function(tilewidth = 10000, verbose = TRUE) {
  say <- function(...) if (verbose) message(...)

  say("Loading mutations")
  gr_tumor <- read_caveman_vcfs(verbose = verbose)
  blacklist <- rtracklayer::import(PATH_BLACKLIST)
  hits <- GenomicRanges::findOverlaps(gr_tumor, blacklist)
  if (length(hits)) gr_tumor <- gr_tumor[-S4Vectors::queryHits(hits)]
  gr_tumor <- gr_tumor[GenomicRanges::seqnames(gr_tumor) != "chrY"]
  say("  ", format(length(gr_tumor), big.mark = ","), " after the blacklist")

  # A mutation whose context could not be resolved - at a contig edge, or in a
  # run of N - has no channel and cannot be modelled.
  unresolved <- is.na(gr_tumor$channel)
  if (any(unresolved)) {
    say("  dropping ", sum(unresolved), " mutations with an unresolvable ",
        "trinucleotide context")
    gr_tumor <- gr_tumor[!unresolved]
  }

  say("Loading copy number")
  gr_copy <- read_ascat_segments(verbose = verbose)

  say("Building CopyTrack at ", tilewidth, " bp")
  gr_CopyTrack <- build_CopyTrack(gr_tumor, gr_copy, tilewidth = tilewidth)

  say("Building SignalTrack at ", tilewidth, " bp")
  gr_SignalTrack <- build_SignalTrack(tilewidth = tilewidth, verbose = verbose)

  say("Standardising covariates")
  covariates <- setdiff(names(GenomicRanges::mcols(gr_SignalTrack)), "bin_weight")
  for (v in covariates) {
    GenomicRanges::mcols(gr_SignalTrack)[[v]] <-
      as.numeric(standardize_covariates(GenomicRanges::mcols(gr_SignalTrack)[[v]]))
  }

  say("Attaching covariates to mutations")
  gr_Mutations <- merge_with_tumor(gr_tumor, gr_SignalTrack, verbose = verbose)

  MutMatrix <- SignaturePPF::getTotalMutations(gr_Mutations)
  SignalTrack <- as.matrix(GenomicRanges::mcols(gr_SignalTrack)[, covariates,
                                                                drop = FALSE])
  CopyTrack <- as.matrix(GenomicRanges::mcols(gr_CopyTrack)[, -1, drop = FALSE])
  CopyTrack <- CopyTrack * gr_CopyTrack$bin_weight
  CopyTrack <- CopyTrack[, colnames(MutMatrix), drop = FALSE]

  stopifnot(nrow(SignalTrack) == nrow(CopyTrack))

  list(gr_Mutations = gr_Mutations,
       SignalTrack = SignalTrack,
       CopyTrack = CopyTrack,
       gr_CopyTrack = gr_CopyTrack,
       gr_SignalTrack = gr_SignalTrack)
}


#' Build one binned cohort object from the raw ICGC tracks
#'
#' @param tilewidth Bin width in bases.
#' @return A list with `gr_Mutations`, `SignalTrack`, `CopyTrack` and the two
#'   bin `GRanges`, ready for [SignaturePPF::SignaturePPF()].
build_icgc_dataset <- function(tilewidth = 2000, verbose = TRUE) {
  say <- function(...) if (verbose) message(...)

  say("Loading mutations")
  gr_tumor <- readRDS(PATH_ICGC_SNV)
  # Blacklisted mutations go now rather than through the bin weights: a mutation
  # inside a blacklisted region has no usable exposure behind it at all.
  blacklist <- rtracklayer::import(PATH_BLACKLIST)
  hits <- GenomicRanges::findOverlaps(gr_tumor, blacklist)
  if (length(hits)) gr_tumor <- gr_tumor[-S4Vectors::queryHits(hits)]
  gr_tumor <- gr_tumor[GenomicRanges::seqnames(gr_tumor) != "chrY"]
  say("  ", format(length(gr_tumor), big.mark = ","), " mutations, ",
      length(unique(gr_tumor$sample)), " samples")

  say("Loading copy number")
  df_copy <- readr::read_tsv(PATH_ICGC_CN, show_col_types = FALSE)
  gr_copy <- GenomicRanges::GRanges(
    seqnames = paste0("chr", df_copy$chr),
    ranges = IRanges::IRanges(start = df_copy$start, end = df_copy$end),
    strand = "*", sample = df_copy$sampleID, score = df_copy$value)

  say("Building CopyTrack at ", tilewidth, " bp")
  gr_CopyTrack <- build_CopyTrack(gr_tumor, gr_copy, tilewidth = tilewidth)

  say("Building SignalTrack at ", tilewidth, " bp")
  gr_SignalTrack <- build_SignalTrack(tilewidth = tilewidth, verbose = verbose)

  say("Standardising covariates")
  covariates <- setdiff(names(GenomicRanges::mcols(gr_SignalTrack)), "bin_weight")
  for (v in covariates) {
    GenomicRanges::mcols(gr_SignalTrack)[[v]] <-
      as.numeric(standardize_covariates(GenomicRanges::mcols(gr_SignalTrack)[[v]]))
  }

  say("Attaching covariates to mutations")
  gr_Mutations <- merge_with_tumor(gr_tumor, gr_SignalTrack, verbose = verbose)

  MutMatrix <- SignaturePPF::getTotalMutations(gr_Mutations)
  SignalTrack <- as.matrix(GenomicRanges::mcols(gr_SignalTrack)[, covariates,
                                                                drop = FALSE])
  # Copy number times usable sequence: this is the exposure per bin, which is
  # what the intensity integrates against.
  CopyTrack <- as.matrix(GenomicRanges::mcols(gr_CopyTrack)[, -1, drop = FALSE])
  CopyTrack <- CopyTrack * gr_CopyTrack$bin_weight
  CopyTrack <- CopyTrack[, colnames(MutMatrix), drop = FALSE]

  stopifnot(nrow(SignalTrack) == nrow(CopyTrack),
            identical(colnames(SignalTrack),
                      intersect(names(GenomicRanges::mcols(gr_Mutations)),
                                colnames(SignalTrack))))

  list(gr_Mutations = gr_Mutations,
       SignalTrack = SignalTrack,
       CopyTrack = CopyTrack,
       gr_CopyTrack = gr_CopyTrack,
       gr_SignalTrack = gr_SignalTrack)
}
