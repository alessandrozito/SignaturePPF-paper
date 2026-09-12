## Produces: no figure. Helper file, sourced by the analysis scripts.

################################################################################
# Helpers for the comparison against TensorSignatures (Vohringer et al.,
# Nat Commun 2021).
#
# WHY THE TWO MODELS ARE NOT TRIVIALLY COMPARABLE
# -----------------------------------------------
# Both extend 96-channel NMF with genomic covariates, but parameterise the
# genomic dependence differently:
#
#                     SignaturePPF                  TensorSignatures
#   genomic effect    continuous log-linear         discrete per-state amplitudes
#                     exp(beta_k' x(t))             on a categorical state tensor
#   resolution        continuous position           fixed categorical states
#   covariates        arbitrary continuous tracks   binarised / discretised states
#   strand asymmetry  not in the model              explicit (transcription and
#                                                   replication strand)
#   copy number       explicit multiplier c_j(t)    not modelled
#   inference         MAP / MCMC, compressive       TensorFlow MLE, rank by AIC
#                     prior selects K
#
# The comparison is therefore run on COMMON GROUND: the ChromHMM 15-state
# annotation. The bins ARE the ChromHMM segments, so the state assignment is
# exact for both methods, and PPF is given the states as one-hot covariates with
# one state dropped as reference. Then:
#
#   * beta_k,state from PPF IS the log enrichment relative to the reference
#     state, and TensorSignatures' state amplitude is the same quantity - no
#     rescaling is needed, unlike the continuous case.
#   * signature spectra are matched one-to-one by cosine similarity.
#   * regional mutation rate is compared as expected vs observed counts in fixed
#     genomic windows.
#
# KNOWN LIMITATION - READ BEFORE REPORTING STRAND RESULTS
# -------------------------------------------------------
# TensorSignatures' headline feature is transcription/replication strand
# asymmetry, which needs to know whether the pyrimidine of each mutation sat on
# the plus or minus strand. `gr_Mutations$channel` has ALREADY been
# pyrimidine-normalised by the preprocessing (it reverse-complements when the
# reference base is A or G), so that orientation is lost and strand() is "*"
# throughout. Recovering it means going back to the raw calls and keeping a flag
# for whether the mutation was reverse-complemented.
#
# Until then both strand axes carry all their mass in the "unknown" index, which
# reduces TS to its genomic-state part - which is exactly the part PPF can also
# express, so the comparison stays apples-to-apples. DO NOT report TS
# strand-asymmetry results from this pipeline.
################################################################################

STD_CHRS <- paste0("chr", c(1:22, "X"))          # no chrY: breast cohorts

# The 15 ChromHMM core states, in ChromHMM's own order.
CHROM_STATES <- c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh",
                  "ZNF.Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC",
                  "ReprPCWk", "Quies")
CHROM_COLS <- c(TssA = "#FF0000", TssAFlnk = "#FF4500", TxFlnk = "#32CD32",
                Tx = "#008000", TxWk = "#006400", EnhG = "#C2E105",
                Enh = "#FFFF00", ZNF.Rpts = "#66CDAA", Het = "#8A91D0",
                TssBiv = "#CD5C5C", BivFlnk = "#E9967A", EnhBiv = "#BDB76B",
                ReprPC = "#808080", ReprPCWk = "#C0C0C0", Quies = "#F0F0F0")


################################################################################
# 1. Usable bases per bin
################################################################################

#' ACGT bases per bin, after removing assembly gaps and blacklisted regions
#'
#' Same convention as the 10 kb loader, so the two datasets share one notion of
#' exposure - except that here the bins are the ChromHMM segments and so have
#' unequal widths.
add_bin_weights <- function(gr, blacklist_file = PATH_BLACKLIST,
                            gaps_file = PATH_GAPS) {
  blacklist <- rtracklayer::import(blacklist_file)
  gap_gr <- rtracklayer::import(gaps_file)

  BinWeight <- GenomicRanges::width(gr)

  # Assembly gaps: count real ACGT instead of the nominal width. Segments
  # spanning centromeres are large, so this is chunked to bound memory.
  ov_gaps <- GenomicRanges::findOverlaps(gr, gap_gr)
  qh <- unique(S4Vectors::queryHits(ov_gaps))
  if (length(qh)) {
    for (chunk in split(qh, ceiling(seq_along(qh) / 2000))) {
      seqs <- BSgenome::getSeq(
        BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr[chunk])
      BinWeight[chunk] <- rowSums(
        Biostrings::alphabetFrequency(seqs)[, c("A", "C", "G", "T"), drop = FALSE])
    }
  }

  # Blacklist: subtract the overlapping width. A bin can meet several blacklist
  # intervals, so the intersections are summed per bin.
  ov <- GenomicRanges::findOverlaps(gr, blacklist)
  if (length(ov)) {
    inter <- GenomicRanges::width(GenomicRanges::pintersect(
      S4Vectors::Pairs(gr[S4Vectors::queryHits(ov)],
                       blacklist[S4Vectors::subjectHits(ov)])))
    black <- tapply(inter, S4Vectors::queryHits(ov), sum)
    idx <- as.integer(names(black))
    BinWeight[idx] <- pmax(BinWeight[idx] - as.numeric(black), 0)
  }
  gr$bin_weight <- BinWeight
  gr
}


################################################################################
# 2. Chromatin-state bins
################################################################################

#' Read a ChromHMM dense BED and normalise its state names
#'
#' "1_TssA" -> "TssA", "8_ZNF/Rpts" -> "ZNF.Rpts".
read_chromhmm <- function(bed_file = PATH_CHROMHMM) {
  if (!file.exists(bed_file)) {
    stop("ChromHMM BED not found: ", bed_file,
         "\n  Set SIGNATUREPPF_CHROMHMM, or see data/README.md.")
  }
  bed <- rtracklayer::import(bed_file)
  bed <- GenomeInfoDb::keepSeqlevels(bed, STD_CHRS, pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(bed) <- STD_CHRS
  GenomeInfoDb::seqlengths(bed) <- GenomeInfoDb::seqlengths(
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[STD_CHRS]
  bed <- GenomicRanges::trim(bed)          # clip anything past the chromosome end

  nm <- gsub("/", ".", sub("^[0-9]+_", "", bed$name))
  if (!all(nm %in% CHROM_STATES)) {
    stop("unexpected state name(s): ",
         paste(setdiff(unique(nm), CHROM_STATES), collapse = ", "))
  }
  bed$state <- factor(nm, levels = CHROM_STATES)
  GenomicRanges::mcols(bed) <- GenomicRanges::mcols(bed)[, "state", drop = FALSE]
  sort(bed)
}


#' Chromatin-state bins with one-hot covariates
#'
#' The bins ARE the ChromHMM segments. The reference state is dropped from the
#' design, so its bins are all-zero and every beta reads directly as the log
#' enrichment relative to that state.
build_chromatin_state_track <- function(bed_file = PATH_CHROMHMM,
                                        reference = "Quies") {
  stopifnot(reference %in% CHROM_STATES)
  bins <- read_chromhmm(bed_file)
  message("segments on chr1-chr22,chrX: ", length(bins))
  bins <- add_bin_weights(bins)
  n0 <- length(bins)
  bins <- bins[bins$bin_weight > 0]
  message("usable bins (bin_weight > 0): ", length(bins),
          "  (dropped ", n0 - length(bins), ")")

  keep <- setdiff(CHROM_STATES, reference)
  X <- matrix(0, length(bins), length(keep), dimnames = list(NULL, keep))
  j <- match(as.character(bins$state), keep)
  ok <- !is.na(j)
  X[cbind(which(ok), j[ok])] <- 1

  GenomicRanges::mcols(bins) <- cbind(
    S4Vectors::DataFrame(bin_weight = bins$bin_weight, state = bins$state),
    methods::as(X, "DataFrame"))
  attr(bins, "reference") <- reference
  bins
}


#' Copy number on the (variable-width) chromatin bins
build_copytrack_states <- function(gr_bins, gr_copy, samples) {
  CN <- matrix(0, nrow = length(gr_bins), ncol = length(samples),
               dimnames = list(NULL, samples))
  gr_copy$score[is.na(gr_copy$score)] <- 2
  gr_copy <- GenomeInfoDb::keepSeqlevels(gr_copy, STD_CHRS, pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(gr_copy) <- STD_CHRS
  GenomeInfoDb::seqlengths(gr_copy) <- GenomeInfoDb::seqlengths(
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[STD_CHRS]

  for (s in seq_along(samples)) {
    gc <- gr_copy[gr_copy$sample == samples[s]]
    if (length(gc) == 0) { CN[, s] <- 1; next }            # no calls -> diploid
    cov_s <- GenomicRanges::coverage(gc, weight = "score")
    sc <- GenomicRanges::binnedAverage(gr_bins, cov_s, varname = "score")$score
    sc[sc < 0.1] <- 0.1                                    # floor, as elsewhere
    CN[, s] <- sc / 2                                      # relative to diploid
  }
  CN
}


#' Assemble a chromatin-state dataset in SignaturePPF form
#'
#' @return A list that [SignaturePPF()] accepts directly, plus the bookkeeping
#'   the TensorSignatures export needs (`bin_of_mut`, `state_of_bin`).
build_chromatin_dataset <- function(gr_tumor, gr_copy, reference = "Quies",
                                    bed_file = PATH_CHROMHMM) {
  gr_states <- build_chromatin_state_track(bed_file, reference)

  gr_tumor <- GenomeInfoDb::keepSeqlevels(gr_tumor, STD_CHRS, pruning.mode = "coarse")
  ov <- GenomicRanges::findOverlaps(gr_tumor, gr_states)
  gr_Mutations <- gr_tumor[S4Vectors::queryHits(ov)]
  bin_of_mut <- S4Vectors::subjectHits(ov)
  samples <- sort(unique(as.character(gr_Mutations$sample)))
  message("mutations kept: ", length(gr_Mutations), " in ", length(samples),
          " samples")

  message("building copy number for ", length(samples), " samples ...")
  CN <- build_copytrack_states(gr_states, gr_copy, samples)

  cov_cols <- setdiff(CHROM_STATES, reference)
  Xbins <- as.matrix(as.data.frame(
    GenomicRanges::mcols(gr_states))[, cov_cols, drop = FALSE])

  # The model reads each mutation's covariates from mcols(gr_Mutations); keep
  # only sample, channel and the design, in the design's own column order.
  meta <- GenomicRanges::mcols(gr_Mutations)[, c("sample", "channel"), drop = FALSE]
  GenomicRanges::mcols(gr_Mutations) <- cbind(
    meta, methods::as(Xbins[bin_of_mut, , drop = FALSE], "DataFrame"))

  # Exposure on the model scale: usable bases x copy number / 2.
  CopyTrack <- CN * gr_states$bin_weight

  list(gr_Mutations = gr_Mutations,
       SignalTrack = Xbins,
       CopyTrack = CopyTrack[, samples, drop = FALSE],
       gr_SignalTrack = gr_states,
       gr_CopyTrack = gr_states,          # states and CN share the same bins
       CopyNumber = CN,
       bin_of_mut = bin_of_mut,
       state_of_bin = gr_states$state,
       reference = reference)
}


################################################################################
# 3. Export the SAME data as a TensorSignatures tensor
#
#    snv : (3, 3, 15, 96, n_samples), all mass in the unknown-strand cell
#    N   : per (state, sample) exposure, broadcast over the 96 channels
#
#    Layout verified against tensorsignatures 0.5.0: TensorSignature.__init__
#    reads samples = snv.shape[-1], p = snv.shape[-2], and the genomic-state
#    dimensions as snv.shape[2:-2]. The state axis comes BEFORE the 96 channels.
#
#    A plain, self-describing set of TSVs is written rather than an HDF5 guessed
#    against one version of the package; python/run_tensorsignatures.py assembles
#    the final input and re-asserts the layout on its side.
################################################################################
export_ts_chromatin <- function(dataChrom, out_dir = DIR_TENSORSIG,
                                tag = "icgc_chromatin") {
  d <- file.path(out_dir, tag)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  cnt <- data.frame(
    tx_strand = 3L, rep_strand = 3L,
    state = as.integer(dataChrom$state_of_bin)[dataChrom$bin_of_mut],
    channel = as.character(dataChrom$gr_Mutations$channel),
    sample = as.character(dataChrom$gr_Mutations$sample),
    stringsAsFactors = FALSE)
  cnt <- dplyr::count(cnt, .data$tx_strand, .data$rep_strand, .data$state,
                      .data$channel, .data$sample, name = "count")

  # Exposure per (state, sample), from the SAME CopyTrack the PPF model uses, so
  # neither method is given a different notion of how much genome each state has.
  expo <- rowsum(dataChrom$CopyTrack, as.integer(dataChrom$state_of_bin))
  expo_df <- as.data.frame(as.table(as.matrix(expo)), stringsAsFactors = FALSE)
  names(expo_df) <- c("state", "sample", "exposure")
  expo_df$state <- as.integer(expo_df$state)

  readr::write_tsv(cnt, file.path(d, "snv_counts_long.tsv.gz"))
  readr::write_tsv(tibble::tibble(state = seq_along(CHROM_STATES),
                                  name = CHROM_STATES,
                                  index = seq_along(CHROM_STATES)),
                   file.path(d, "state_key.tsv"))
  readr::write_tsv(expo_df, file.path(d, "state_exposure.tsv"))
  writeLines(c(
    sprintf("dims: (tx=3, rep=3, state=%d, channel=96, sample=%d)",
            length(CHROM_STATES), length(unique(cnt$sample))),
    "state axis: the 15 ChromHMM core states (see state_key.tsv);",
    "  the bins ARE the ChromHMM segments, so the assignment is exact",
    "strand axes: uninformative (all mass at index 3) by design - see the",
    "  limitation note at the top of R/TensorSignatures_functions.R",
    "state_exposure.tsv: bin_weight * CN/2 summed per (state, sample)",
    sprintf("total counts: %d", sum(cnt$count))),
    file.path(d, "README.txt"))
  message("written to ", d)
  invisible(d)
}


################################################################################
# 4. Read a TensorSignatures fit
################################################################################
read_ts_fit <- function(dir) {
  f_sig <- file.path(dir, "ts_signatures.tsv")
  f_amp <- file.path(dir, "ts_state_amplitudes.tsv")
  f_key <- file.path(dir, "state_key.tsv")
  if (!file.exists(f_sig)) {
    warning("no TensorSignatures output in ", dir,
            " - run the sweep first (R/Comparison_TensorSignatures.R)")
    return(NULL)
  }
  sig <- as.data.frame(readr::read_tsv(f_sig, show_col_types = FALSE))
  rownames(sig) <- sig[[1]]          # first column is the channel name
  sig[[1]] <- NULL
  list(signatures = as.matrix(sig),
       amplitudes = if (file.exists(f_amp))
         readr::read_tsv(f_amp, show_col_types = FALSE) else NULL,
       state_key = if (file.exists(f_key))
         readr::read_tsv(f_key, show_col_types = FALSE) else NULL,
       dir = dir)
}


#' Rank selection: collect the fit summaries a sweep produced
ts_sweep_summary <- function(base_dir) {
  dirs <- list.files(base_dir, pattern = "^rank[0-9]+$", full.names = TRUE)
  if (!length(dirs)) {
    warning("no rank* directories in ", base_dir)
    return(NULL)
  }
  out <- lapply(dirs, function(d) {
    f <- file.path(d, "ts_fit_summary.tsv")
    if (!file.exists(f)) return(NULL)
    s <- readr::read_tsv(f, show_col_types = FALSE)
    s$dir <- d
    s
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  out[order(out$rank), ]
}


################################################################################
# 5. Comparison
################################################################################

#' PPF chromatin-state effects
#'
#' With one-hot covariates, beta_k,state IS the log enrichment relative to the
#' reference state, so no rescaling is needed (unlike the continuous case, where
#' the effect has to be multiplied by the covariate contrast between states).
ppf_chromatin_effects <- function(fit) {
  df <- as.data.frame(as.table(fit$Betas), stringsAsFactors = FALSE)
  names(df) <- c("state", "signature", "beta")
  df$state <- as.character(df$state)
  df$signature <- as.character(df$signature)
  df
}


#' PPF effect of a CONTINUOUS covariate, on the scale of a two-state contrast
#'
#' For a median split of covariate l, the model implies a log intensity ratio
#' between the states of beta_kl * (E[x_l | high] - E[x_l | low]). That is the
#' quantity TensorSignatures reports as a state amplitude, so this puts the two
#' on the same scale. Only needed for the continuous-covariate variant of the
#' comparison; the chromatin-state version is already on that scale.
ppf_state_logratio <- function(fit, state_means) {
  B <- fit$Betas
  vars <- intersect(rownames(B), names(state_means))
  do.call(rbind, lapply(vars, function(v) {
    m <- state_means[[v]]
    dx <- as.numeric(m[length(m)] - m[1])
    data.frame(covariate = v, signature = colnames(B),
               beta = as.numeric(B[v, ]), delta_x = dx,
               ppf_logratio = as.numeric(B[v, ]) * dx,
               stringsAsFactors = FALSE)
  }))
}


#' One-to-one matching between two signature sets
#'
#' Maximum total cosine similarity by Hungarian assignment. This replaces
#' per-signature argmax matching, which is many-to-one: with argmax two TS
#' signatures can both claim the same PPF signature - duplicating its rows -
#' while other PPF signatures are never claimed and vanish in the join. The
#' Hungarian solution gives min(K_ts, K_ppf) unique pairs; anything left over is
#' returned in the attributes rather than disappearing.
hungarian_match_signatures <- function(S_ts, S_ppf) {
  S_ppf <- S_ppf[rownames(S_ts), , drop = FALSE]
  cosine <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))
  cs <- outer(seq_len(ncol(S_ts)), seq_len(ncol(S_ppf)),
              Vectorize(function(a, b) cosine(S_ts[, a], S_ppf[, b])))
  pairs <- RcppHungarian::HungarianSolver(1 - cs)$pairs   # maximise total cosine
  pairs <- pairs[pairs[, 1] > 0 & pairs[, 2] > 0, , drop = FALSE]
  out <- tibble::tibble(signature = colnames(S_ts)[pairs[, 1]],   # TS label
                        ppf_sig = colnames(S_ppf)[pairs[, 2]],
                        cosine = cs[pairs])
  attr(out, "unmatched_ts") <- setdiff(colnames(S_ts), out$signature)
  attr(out, "unmatched_ppf") <- setdiff(colnames(S_ppf), out$ppf_sig)
  out
}


################################################################################
# 5b. Labels
#
#    Neither method labels its signatures meaningfully. PPF numbers them in
#    whatever order the initialisation happened to put them; TensorSignatures in
#    whatever order TensorFlow happened to converge to. Two panels side by side
#    are only readable if "01" means the same thing in both, so both are
#    renumbered: PPF by decreasing relevance weight with relabel_by_mu() in
#    Utils_functions.R, and TensorSignatures by which PPF signature it matches,
#    below.
################################################################################

#' Match TensorSignatures to a mu-ordered PPF fit and renumber it to agree
#'
#' `TS01` becomes whichever TS signature the Hungarian assignment pairs with
#' `SigN01`, `TS02` the partner of `SigN02`, and so on - so the two methods can
#' be drawn side by side with the panels lining up.
#'
#' Dead PPF signatures are dropped BEFORE the matching, not after. The assignment
#' is one-to-one, so a signature the compressive prior has switched off can win a
#' TS signature and displace a live one; filtering afterwards would leave that TS
#' signature unpaired rather than paired correctly. Its coefficients are draws
#' from the prior in any case, so comparing them against a TS amplitude would
#' measure the prior.
#'
#' @param ts Output of [read_ts_fit()].
#' @param fit A fit already renumbered by [relabel_by_mu()] - the live signatures
#'   have to be a prefix for the TS numbering to come out contiguous.
#' @param mu_min Relevance weight below which a PPF signature counts as dead.
#' @return A list of `ts` (renamed, and with its signature columns in the new
#'   order), `match` (the matching table, carrying both TS labels), and `live`.
match_and_relabel_ts <- function(ts, fit, mu_min = 0.05, prefix = "TS") {
  live <- colnames(fit$Signatures)[fit$Mu >= mu_min]
  if (!length(live)) {
    stop("every PPF signature has mu < ", mu_min, ": nothing left to compare.",
         call. = FALSE)
  }
  if (length(live) < ncol(fit$Signatures)) {
    message("dropping ", ncol(fit$Signatures) - length(live), " of ",
            ncol(fit$Signatures), " PPF signature(s) with mu < ", mu_min, ": ",
            paste(setdiff(colnames(fit$Signatures), live), collapse = ", "))
  }

  map <- hungarian_match_signatures(ts$signatures,
                                    fit$Signatures[, live, drop = FALSE])
  unmatched_ppf <- attr(map, "unmatched_ppf")

  # Index of each TS signature = index of its PPF partner. Anything TS found that
  # PPF did not is numbered after the live block, in TensorSignatures' own order,
  # so it is visibly an extra rather than silently occupying a matched slot.
  idx   <- stats::setNames(match(map$ppf_sig, live), map$signature)
  extra <- setdiff(colnames(ts$signatures), names(idx))
  if (length(extra)) idx[extra] <- length(live) + seq_along(extra)
  lab <- stats::setNames(sprintf("%s%02d", prefix, idx), names(idx))

  o <- order(idx[colnames(ts$signatures)])
  old_cols <- colnames(ts$signatures)[o]
  ts$signatures <- ts$signatures[, o, drop = FALSE]
  colnames(ts$signatures) <- unname(lab[old_cols])

  # `signature` stays the join key against the amplitude table, so it has to hold
  # the NEW label on both sides; the original is kept as `ts_original` because it
  # is the only way back to the raw TensorSignatures output on disk.
  if (!is.null(ts$amplitudes)) {
    ts$amplitudes$signature <- unname(lab[ts$amplitudes$signature])
  }
  map$ts_original <- map$signature
  map$signature   <- unname(lab[map$ts_original])
  map <- map[order(match(map$ppf_sig, live)),
             c("signature", "ppf_sig", "cosine", "ts_original")]

  attr(map, "unmatched_ts")  <- unname(lab[extra])
  attr(map, "unmatched_ppf") <- unmatched_ppf
  attr(map, "dead_ppf")      <- setdiff(colnames(fit$Signatures), live)

  list(ts = ts, match = map, live = live)
}


#' Put both methods' chromatin-state effects on one scale and join them
#'
#' TS amplitudes are relative to ITS state 1, so they are re-referenced to the
#' PPF reference state before joining.
#'
#' @param ts Output of [read_ts_fit()], already renumbered by
#'   [match_and_relabel_ts()].
#' @param map The matching table from [match_and_relabel_ts()]. Passed in rather
#'   than recomputed here: the matching depends on which PPF signatures are
#'   considered live, so a second run of it inside this function could disagree
#'   with the one the labels were built from, and the labels would then be a
#'   quiet lie.
compare_chromatin_effects <- function(fit, ts, map, reference = "Quies") {
  if (is.null(ts)) return(NULL)

  # `ppf_sig`, deliberately not `signature`: the amplitude table already has a
  # `signature` column (the TS label) and a clash silently yields
  # signature.x / signature.y.
  live <- map$ppf_sig
  fit$Signatures <- fit$Signatures[, live, drop = FALSE]
  fit$Betas <- fit$Betas[, live, drop = FALSE]
  fit$Mu <- fit$Mu[live]
  amp <- ts$amplitudes

  # The Python side already writes a readable state `name`; only join the key
  # when it is missing, or the two `name` columns collide.
  if (!"name" %in% names(amp)) {
    amp <- dplyr::left_join(amp,
                            dplyr::select(ts$state_key, .data$state, .data$name),
                            by = "state")
  }

  # Rename BEFORE transmute: inside transmute() expressions are evaluated in
  # order, so `ts_signature = signature` written after `signature = ppf_sig`
  # would pick up the NEW column and duplicate ppf_sig into both.
  ts_eff <- amp |>
    dplyr::left_join(map, by = "signature") |>
    dplyr::rename(ts_signature = "signature") |>
    dplyr::transmute(state = as.character(.data$name),
                     signature = .data$ppf_sig,
                     ts_signature = .data$ts_signature,
                     cosine = .data$cosine,
                     ts_logratio = .data$log_ratio_vs_state1)

  ts_eff$pair_label <- sprintf("%s ~ %s  (%.2f)", ts_eff$signature,
                               ts_eff$ts_signature, ts_eff$cosine)

  # Re-reference the TS amplitudes from ITS state 1 to the PPF reference state.
  ref <- ts_eff[ts_eff$state == reference, c("ts_signature", "ts_logratio")]
  names(ref)[2] <- "r"
  ts_eff <- dplyr::left_join(ts_eff, ref, by = "ts_signature")
  ts_eff$ts_logratio <- ts_eff$ts_logratio - ts_eff$r
  ts_eff$r <- NULL

  res <- dplyr::inner_join(ppf_chromatin_effects(fit), ts_eff,
                           by = c("state", "signature"))
  # Panels in label order, i.e. by decreasing mu, NOT by decreasing cosine: the
  # whole point of the renumbering is that panel i is the same process wherever
  # it is drawn, so the order has to be a property of the signature rather than
  # of how well it happened to match.
  lev <- unique(res$pair_label[order(res$signature)])
  res$pair_label <- factor(res$pair_label, levels = lev)
  attr(res, "unmatched_ts") <- attr(map, "unmatched_ts")
  attr(res, "unmatched_ppf") <- attr(map, "unmatched_ppf")
  res
}


################################################################################
# 6. Regional mutation rate
#
#    PPF gives an intensity per bin directly.
#
#    TensorSignatures only predicts a total per (state, sample) - it has no
#    notion of position WITHIN a state - so its per-bin prediction is that total
#    spread over the segments of the state in proportion to each segment's
#    exposure. This is not a handicap imposed by us: it is the finest resolution
#    the model can express, and it is exactly the difference the comparison is
#    meant to expose.
#
#    Both predictions and the observed counts are then aggregated into fixed
#    genomic windows so the comparison is on a common, interpretable scale rather
#    than on the unequal ChromHMM segments.
################################################################################
compare_mutation_rate <- function(dataChrom, fit, ts_dir, window = 1e6) {
  f <- file.path(ts_dir, "ts_predicted_state_sample.tsv")
  if (!file.exists(f)) {
    stop("no ts_predicted_state_sample.tsv in ", ts_dir,
         " - refit with the current python/run_tensorsignatures.py")
  }
  samples <- colnames(dataChrom$CopyTrack)

  n_bins <- length(dataChrom$gr_SignalTrack)
  obs <- unclass(table(
    factor(dataChrom$bin_of_mut, levels = seq_len(n_bins)),
    factor(as.character(dataChrom$gr_Mutations$sample), levels = samples)))
  storage.mode(obs) <- "double"

  # PPF intensity per bin. reconstruct_lambda() takes the BASELINE from the fit,
  # which under the activity prior is $Baseline and NOT $Thetas.
  Lam_ppf <- reconstruct_lambda(fit, dataChrom$SignalTrack,
                                dataChrom$CopyTrack[, samples, drop = FALSE])

  # TS per-state total, spread over that state's bins by exposure share.
  pred <- readr::read_tsv(f, show_col_types = FALSE)
  P <- matrix(0, length(CHROM_STATES), length(samples),
              dimnames = list(NULL, samples))
  P[cbind(pred$state, match(pred$sample, samples))] <- pred$predicted
  st <- as.integer(dataChrom$state_of_bin)
  tot_expo <- rowsum(dataChrom$CopyTrack[, samples, drop = FALSE], st)
  share <- dataChrom$CopyTrack[, samples, drop = FALSE] /
    tot_expo[st, , drop = FALSE]
  Lam_ts <- P[st, , drop = FALSE] * share

  gr <- dataChrom$gr_SignalTrack
  win <- paste0(as.character(GenomicRanges::seqnames(gr)), ":",
                floor((GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2 / window))
  win <- factor(win, levels = unique(win))          # genomic order preserved
  O <- rowsum(obs, win); A <- rowsum(Lam_ppf, win); B <- rowsum(Lam_ts, win)

  data.frame(index = seq_len(nrow(O)),
             window = rownames(O),
             chrom = sub(":.*", "", rownames(O)),
             observed = rowSums(O), ppf = rowSums(A), tensorsig = rowSums(B),
             row.names = NULL)
}


#' Headline numbers: how well each model predicts regional burden
score_mutation_rate <- function(rate) {
  data.frame(
    model = c("SignaturePPF", "TensorSignatures"),
    rmse = c(sqrt(mean((rate$ppf - rate$observed)^2)),
             sqrt(mean((rate$tensorsig - rate$observed)^2))),
    cor = c(cor(rate$ppf, rate$observed), cor(rate$tensorsig, rate$observed)),
    total_predicted = c(sum(rate$ppf), sum(rate$tensorsig)),
    total_observed = sum(rate$observed),
    row.names = NULL)
}


################################################################################
# 7. Figures
################################################################################

plot_chromatin_betas <- function(fit, reference = "Quies") {
  df <- ppf_chromatin_effects(fit)
  df$state <- factor(df$state, levels = setdiff(CHROM_STATES, reference))
  df$signature <- factor(df$signature, levels = colnames(fit$Betas))
  ggplot2::ggplot(df, ggplot2::aes(.data$state, .data$beta, fill = .data$state)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ signature, scales = "free_y") +
    ggplot2::scale_fill_manual(values = CHROM_COLS, guide = "none") +
    ggplot2::labs(x = NULL,
                  y = bquote(beta ~ "(log enrichment vs" ~ .(reference) * ")"),
                  title = "Chromatin-state effects estimated by SignaturePPF") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = 7))
}


#' PPF beta against TensorSignatures amplitude
#'
#' The same points, cut two ways. Both cuts are needed because a pooled
#' correlation can look reasonable while the structure underneath disagrees:
#'
#' * `by = "signature"` gives one panel per matched pair, each showing how one
#'   signature's profile ACROSS the states compares. This asks whether the two
#'   methods agree on where a given process is enriched.
#' * `by = "state"` gives one panel per chromatin state, each showing how the
#'   signatures compare WITHIN that state. This asks whether the two methods
#'   agree on what a given piece of chromatin does, and it is the cut that
#'   exposes a state on which the methods disagree for every signature at once -
#'   which in the by-signature view is spread thin across every panel and easy
#'   to miss.
#' * `by = "none"` pools everything into one panel.
#'
#' The colour follows the cut: faceting by signature colours by state and vice
#' versa, so in either view the point's other coordinate is still readable.
#'
#' @param cmp Output of [compare_chromatin_effects()].
#' @param by What each panel holds fixed.
#' @param ncol Panels per row.
#' @param free_scales Free axes per panel. The states span a far wider range of
#'   effects than the signatures do, so fixed scales leave most panels empty.
plot_chromatin_effects <- function(cmp, by = c("signature", "state", "none"),
                                   ncol = 4, free_scales = TRUE) {
  by <- match.arg(by)
  cmp$state <- factor(cmp$state, levels = CHROM_STATES)
  sigs <- sort(unique(as.character(cmp$signature)))
  cmp$signature <- factor(cmp$signature, levels = sigs)

  facet_var <- switch(
    by,
    none = NULL,
    state = "state",
    signature = if ("pair_label" %in% names(cmp)) "pair_label" else "signature")

  p <- ggplot2::ggplot(cmp, ggplot2::aes(.data$ts_logratio, .data$beta)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70")

  p <- if (identical(by, "state")) {
    p + ggplot2::geom_point(ggplot2::aes(colour = .data$signature), size = 1.8) +
      ggplot2::scale_colour_manual(values = sig_palette(sigs), drop = FALSE,
                                   name = "Signature")
  } else {
    p + ggplot2::geom_point(ggplot2::aes(colour = .data$state), size = 1.8) +
      ggplot2::scale_colour_manual(values = CHROM_COLS, drop = FALSE,
                                   name = "State")
  }

  p <- p +
    ggplot2::labs(x = "TensorSignatures log amplitude",
                  y = expression(SignaturePPF ~ beta)) +
    ggplot2::theme_bw()

  if (is.null(facet_var)) return(p)

  p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), ncol = ncol,
                          scales = if (free_scales) "free" else "fixed")
}


#' Observed vs predicted burden along the genome, with chromosome bands
plot_mutation_rate <- function(rate,
                               cols = c(SignaturePPF = "#CD2626",
                                        TensorSignatures = "#000D8B")) {
  chr <- unique(rate$chrom)
  band <- stats::setNames(seq_along(chr) %% 2, chr)
  rects <- do.call(rbind, lapply(split(rate, rate$chrom), function(z) {
    data.frame(chrom = z$chrom[1], xmin = min(z$index), xmax = max(z$index) + 1)
  }))
  rects$fill <- factor(band[rects$chrom])

  long <- rbind(
    data.frame(index = rate$index, model = "SignaturePPF", pred = rate$ppf),
    data.frame(index = rate$index, model = "TensorSignatures", pred = rate$tensorsig))

  ggplot2::ggplot() +
    ggplot2::geom_rect(data = rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = 0, ymax = Inf, fill = .data$fill),
                       alpha = 0.10, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("0" = "#4682B4", "1" = "antiquewhite")) +
    ggplot2::geom_point(data = rate,
                        ggplot2::aes(.data$index, .data$observed),
                        size = 0.45, colour = "grey35") +
    ggplot2::geom_line(data = long,
                       ggplot2::aes(.data$index, .data$pred, colour = .data$model),
                       linewidth = 0.6, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = "Genomic window", y = "Mutations",
                  title = "Regional mutation rate: observed (grey) vs predicted") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
