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
