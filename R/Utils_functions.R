## Adapters between the preprocessed cohort objects and SignaturePPF.

#' Coerce a preprocessed cohort object into the form SignaturePPF expects
#'
#' The objects saved by the old loaders carry more than the model needs -
#' `gr_CopyTrack`, `gr_SignalTrack`, and a `tumor` column in the mutation
#' metadata. SignaturePPF_validate() ignores the extra list elements, but it
#' reads the covariates out of `mcols(gr_Mutations)` by taking every numeric
#' column that is not `sample` or `channel`, so any stray numeric column would be
#' silently fitted as a genomic covariate.
#'
#' This restricts the metadata to exactly `sample`, `channel` and the columns of
#' `SignalTrack`, IN THE ORDER SignalTrack has them, which is the alignment the
#' model depends on and cannot check for itself.
#'
#' @param data A list with gr_Mutations, SignalTrack and CopyTrack.
#' @param covariates Covariates to keep. Defaults to all of `SignalTrack`.
#' @param drop_samples Samples to exclude from both the mutations and CopyTrack.
#' @return A list with gr_Mutations, SignalTrack, CopyTrack, plus the untouched
#'   `gr_SignalTrack` / `gr_CopyTrack` when present (the TensorSignatures step
#'   needs the bin ranges).
as_ppf_data <- function(data, covariates = NULL, drop_samples = character(0)) {
  stopifnot(all(c("gr_Mutations", "SignalTrack", "CopyTrack") %in% names(data)))

  gr <- data$gr_Mutations
  SignalTrack <- as.matrix(data$SignalTrack)
  CopyTrack <- as.matrix(data$CopyTrack)

  if (is.null(covariates)) covariates <- colnames(SignalTrack)
  missing_cov <- setdiff(covariates, colnames(SignalTrack))
  if (length(missing_cov)) {
    stop("SignalTrack has no column(s): ", paste(missing_cov, collapse = ", "))
  }
  missing_mc <- setdiff(covariates, colnames(GenomicRanges::mcols(gr)))
  if (length(missing_mc)) {
    stop("mcols(gr_Mutations) has no column(s): ",
         paste(missing_mc, collapse = ", "))
  }

  if (length(drop_samples)) {
    keep_mut <- !(as.character(gr$sample) %in% drop_samples)
    gr <- gr[keep_mut]
    # droplevels matters: an unused factor level would still appear as a column
    # of the mutation matrix and then demand a CopyTrack column that is gone.
    gr$sample <- droplevels(as.factor(gr$sample))
    CopyTrack <- CopyTrack[, !(colnames(CopyTrack) %in% drop_samples), drop = FALSE]
  }

  GenomicRanges::mcols(gr) <-
    GenomicRanges::mcols(gr)[, c("sample", "channel", covariates), drop = FALSE]

  out <- list(gr_Mutations = gr,
              SignalTrack = SignalTrack[, covariates, drop = FALSE],
              CopyTrack = CopyTrack)
  for (nm in c("gr_SignalTrack", "gr_CopyTrack")) {
    if (!is.null(data[[nm]])) out[[nm]] <- data[[nm]]
  }
  out
}


#' Load a cohort and put it in model form in one step
#'
#' @param path Path to the preprocessed `.rds.gzip`.
#' @param ... Passed to [as_ppf_data()].
load_cohort <- function(path, ...) {
  if (!file.exists(path)) {
    stop("cohort file not found: ", path,
         "\n  Set SIGNATUREPPF_DATA, or see data/README.md.")
  }
  message("loading ", basename(path))
  as_ppf_data(readRDS(path), ...)
}


#' Per-bin, per-sample expected mutation count under a fitted model
#'
#' \deqn{\Lambda_{bj} = c_j(b) \sum_k \phi_{kj} e^{\beta_k' x(b)}}
#'
#' The old package computed this in C++ (`Reconstruct_Lambda`), but the sum over
#' channels collapses - the signature columns are normalised - so it reduces to
#' one matrix product and does not need compiled code.
#'
#' NOTE the second argument: this needs the BASELINE \eqn{\phi}, not the total
#' activity. Under the activity prior `fit$Thetas` is \eqn{\theta = \phi q_j(\beta_k)}
#' and passing it would inflate every intensity by \eqn{q_j(\beta_k)}.
#'
#' @param fit A `SignaturePPF` fit.
#' @param SignalTrack,CopyTrack The bins the model was fitted on.
#' @return An `nbins x nsamples` matrix of expected counts.
reconstruct_lambda <- function(fit, SignalTrack, CopyTrack) {
  samples <- colnames(CopyTrack)
  Phi <- fit$Baseline[, samples, drop = FALSE]
  E <- exp(pmin(pmax(SignalTrack[, rownames(fit$Betas), drop = FALSE] %*% fit$Betas,
                     -20), 20))
  CopyTrack * (E %*% Phi)
}


#' Bin index of every mutation
#'
#' Needs `gr_SignalTrack`, the bin ranges, which the preprocessed cohort objects
#' carry alongside the covariate matrix. The bins tile the retained genome
#' without overlapping, so every mutation hits at most one.
#'
#' A mutation that hits NO bin is returned as `NA` rather than dropped, so the
#' caller decides. Any subsetting of the form `bin_of_mut %in% bins` excludes
#' them automatically, which is the behaviour wanted here: a mutation in an
#' excluded region belongs to neither the training nor the held-out set.
#'
#' @param data A cohort object with `gr_Mutations` and `gr_SignalTrack`.
#' @return An integer vector, one bin index per mutation.
bin_of_mutation <- function(data) {
  if (is.null(data$gr_SignalTrack)) {
    stop("`data` has no `gr_SignalTrack`, so mutations cannot be mapped onto ",
         "bins. Load the cohort with load_cohort(), which keeps it.",
         call. = FALSE)
  }
  idx <- GenomicRanges::findOverlaps(data$gr_Mutations, data$gr_SignalTrack,
                                     select = "first")
  if (anyNA(idx)) {
    message(sum(is.na(idx)), " of ", length(idx),
            " mutations fall outside every bin and are excluded.")
  }
  idx
}


#' Restrict a cohort to a set of bins and a set of covariates
#'
#' Used to fit on a training subset of the genome, and to fit the sequence of
#' nested covariate sets in the stability analysis.
#'
#' The covariates are NOT re-standardised on the subset. They were standardised
#' over the whole genome, and the values carried at the mutations were
#' standardised with them; rescaling here would put the bins and the mutations on
#' different scales, and would also put the training and held-out bins on
#' different scales, making a coefficient fitted on one meaningless on the other.
#'
#' @param data A cohort object, as returned by [load_cohort()].
#' @param bins Integer indices of the bins to keep. Sorted on the way in: the
#'   rows of `SignalTrack` and `CopyTrack` must stay in genomic order.
#' @param bin_of_mut Bin index of each mutation, from [bin_of_mutation()].
#' @param covariates Covariates to keep. Defaults to all of them.
#' @return A cohort object holding only those bins, and only the mutations
#'   inside them.
subset_bins <- function(data, bins, bin_of_mut, covariates = NULL) {
  if (is.null(covariates)) covariates <- colnames(data$SignalTrack)
  bins <- sort(unique(as.integer(bins)))

  gr <- data$gr_Mutations[which(bin_of_mut %in% bins)]
  GenomicRanges::mcols(gr) <-
    GenomicRanges::mcols(gr)[, c("sample", "channel", covariates), drop = FALSE]
  # A sample left with no mutations would otherwise survive as an all-zero
  # column of the mutation matrix and be given an activity fitted to nothing.
  gr$sample <- droplevels(as.factor(gr$sample))

  list(gr_Mutations = gr,
       SignalTrack = data$SignalTrack[bins, covariates, drop = FALSE],
       CopyTrack = data$CopyTrack[bins, , drop = FALSE],
       gr_SignalTrack = data$gr_SignalTrack[bins])
}


#' Expected counts of a fitted model on an arbitrary set of bins
#'
#' The bins need not be the ones the model was fitted on: with \eqn{\phi} held at
#' its fitted value the intensity extends to any bin whose covariates and copy
#' number are known, which is what makes an out-of-sample prediction possible.
#'
#' @param fit A `SignaturePPF` fit.
#' @param data The cohort the bins are indexed into.
#' @param bins Integer bin indices.
#' @return A `length(bins) x nsamples` matrix of expected counts.
predict_lambda_bins <- function(fit, data, bins) {
  reconstruct_lambda(
    fit,
    data$SignalTrack[bins, , drop = FALSE],
    data$CopyTrack[bins, colnames(fit$Thetas), drop = FALSE])
}


#' Per-patient RMSE of the predicted mutation rate, at a coarser resolution
#'
#' Both prediction and observation are aggregated from the model's bins up to
#' wider regions before comparing. At 10 kb the counts are mostly 0 and 1 and the
#' RMSE is dominated by Poisson noise no model can predict; aggregating to 1 Mb
#' asks the question actually of interest, whether the model gets the REGIONAL
#' rate right.
#'
#' @param Lambda Expected counts on `bins`, from [predict_lambda_bins()].
#' @param obs Observed counts, bins in rows and samples in columns, over the
#'   whole genome.
#' @param bins The bins `Lambda` was computed on.
#' @param region_of_bin Region index of every bin in the genome.
#' @return A named vector of RMSEs, one per patient.
patient_rmse <- function(Lambda, obs, bins, region_of_bin) {
  grp <- region_of_bin[bins]
  keep <- !is.na(grp)
  agg_pred <- rowsum(Lambda[keep, , drop = FALSE], grp[keep])
  agg_obs <- rowsum(obs[bins[keep], colnames(Lambda), drop = FALSE], grp[keep])
  sqrt(colMeans((agg_pred - agg_obs)^2))
}


#' Observed counts per bin and sample
#'
#' @param bin_of_mut Bin index of each mutation.
#' @param sample_of_mut Sample of each mutation.
#' @param n_bins Total number of bins.
#' @return An `n_bins x nsamples` integer matrix.
count_by_bin <- function(bin_of_mut, sample_of_mut, n_bins) {
  tab <- table(factor(bin_of_mut, levels = seq_len(n_bins)),
               as.character(sample_of_mut))
  # as.integer(), not as.matrix(): a table of this size is worth keeping in four
  # bytes an entry rather than eight.
  matrix(as.integer(tab), nrow = n_bins,
         dimnames = list(NULL, colnames(tab)))
}


#' Best cosine similarity of each column against a reference catalogue
#'
#' Replaces `SigPoisProcess::match_to_RefSigs()`, which was not carried over into
#' SignaturePPF.
#'
#' @param sigs A signature matrix, channels in rows.
#' @param ref Reference catalogue. Defaults to COSMIC v3.4 SBS96.
#' @return A data frame with the best-matching reference and its cosine.
match_to_cosmic <- function(sigs, ref = SignaturePPF::COSMIC_v3.4_SBS96_GRCh37) {
  sigs <- as.matrix(sigs)[rownames(ref), , drop = FALSE]
  cs <- crossprod(sigs, ref) /
    outer(sqrt(colSums(sigs^2)), sqrt(colSums(ref^2)))
  best <- max.col(cs, ties.method = "first")
  data.frame(signature = colnames(sigs),
             best_match = colnames(ref)[best],
             cosine = cs[cbind(seq_len(nrow(cs)), best)],
             row.names = NULL, stringsAsFactors = FALSE)
}
## Assigning individual mutations to signatures under a fitted model.

#' The baseline phi_kj of a fit
#'
#' The assignment probability of a mutation to signature k is proportional to
#' \eqn{r_{ik}\,\phi_{kj}\,e^{\beta_k'x(t)}}, so it is the BASELINE that enters,
#' not the activity.
#'
#' SignaturePPF fits the activity prior, where `$Thetas` holds the total activity
#' \eqn{\theta_{kj} = \phi_{kj} q_j(\beta_k)} and `$Baseline` holds \eqn{\phi}.
#' Fits from the predecessor package's original prior have no `$Baseline` at all
#' and their `$Thetas` IS \eqn{\phi}. Both are accepted so that old fits can be
#' re-analysed, because getting this wrong is silent: \eqn{q_j(\beta_k)} varies
#' with k, so it does not cancel in the argmax and the assignment simply changes.
get_baseline <- function(fit) {
  if (!is.null(fit$Baseline)) fit$Baseline else fit$Thetas
}


#' Most probable signature for every mutation
#'
#' @param fit A fitted model.
#' @param data The cohort it was fitted on, or a bare `GRanges` of mutations.
#'   The second form is for comparing SEVERAL fits on ONE fixed set of
#'   mutations: `mcols` need only carry the covariates that fit uses, so a set
#'   of mutations can be pushed through a whole sequence of nested models.
#' @return A character vector, one signature name per mutation.
assign_mutations <- function(fit, data) {
  gr <- if (methods::is(data, "GRanges")) data else data$gr_Mutations
  Phi <- get_baseline(fit)
  X <- as.matrix(GenomicRanges::mcols(gr)[, rownames(fit$Betas), drop = FALSE])

  # as.character() is not decoration: `channel` and `sample` are factors, and
  # indexing a matrix with a factor uses its integer CODES, not its labels. That
  # happens to be right when the fit was built on this exact object and wrong,
  # silently, as soon as the signatures have been pruned or reordered.
  bigProd <- fit$Signatures[as.character(gr$channel), , drop = FALSE] *
    t(Phi[, as.character(gr$sample), drop = FALSE]) *
    exp(pmin(pmax(X %*% fit$Betas, -20), 20))

  # max.col, not apply(., 1, which.max): identical result including ties, but it
  # does not loop in R over several hundred thousand mutations.
  colnames(fit$Signatures)[max.col(bigProd, ties.method = "first")]
}


#' Mutations attributed to each signature, alongside its relevance weight
#'
#' @param fit A fitted model.
#' @param data The cohort it was fitted on.
#' @return A data frame with one row per signature: `best_sig`, the number of
#'   mutations `m` assigned to it, and its `mu`. Signatures that win no mutation
#'   are kept with `m = 0` - they are exactly the ones worth seeing.
df_assign <- function(fit, data) {
  levs <- colnames(fit$Signatures)
  counts <- table(factor(assign_mutations(fit, data), levels = sort(levs)))
  data.frame(best_sig = names(counts), m = as.integer(counts),
             row.names = NULL, stringsAsFactors = FALSE) |>
    merge(data.frame(best_sig = names(fit$Mu), mu = as.numeric(fit$Mu),
                     stringsAsFactors = FALSE),
          by = "best_sig", all.x = TRUE)
}
