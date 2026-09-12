## Helper file, sourced by the analysis scripts.

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


#' Expected count matrix under a fit
#'
#' The `I x J` matrix of expected mutation counts per channel and sample. Under
#' the activity prior this needs no pass over the genome at all:
#' \deqn{M_{ij} = \sum_b \tfrac12 c_j(b) \sum_k r_{ik}\phi_{kj}e^{\beta_k'x_b}
#'             = \sum_k r_{ik}\phi_{kj} q_j(\beta_k) = \sum_k r_{ik}\theta_{kj},}
#' because \eqn{q_j(\beta_k)} is exactly the integral the sum over bins performs.
#' So the whole reconstruction is `Signatures %*% Thetas` - which is the sense in
#' which the activity parametrisation contains an ordinary NMF, and the reason
#' the predecessor's genome-wide `Reconstruct_CountMatrix()` has no counterpart
#' here. `SignalTrack` and `CopyTrack` are accepted and ignored, so call sites
#' ported from that function keep reading the same way.
#'
#' @param fit A `SignaturePPF` fit.
#' @param SignalTrack,CopyTrack Unused; see above.
#' @return An `I x J` matrix of expected counts.
reconstruct_count_matrix <- function(fit, SignalTrack = NULL, CopyTrack = NULL) {
  fit$Signatures %*% fit$Thetas
}


#' The covariate columns of a mutation set
#'
#' `mcols(gr_Mutations)` carries the three identifier columns first and the
#' design after; the model reads each mutation's covariates from there.
mutation_covariates <- function(gr_Mutations) {
  stopifnot(identical(names(GenomicRanges::mcols(gr_Mutations))[1:3],
                      c("tumor", "sample", "channel")))
  as.matrix(GenomicRanges::mcols(gr_Mutations)[, -c(1:3)])
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


#' Donor-level clinical annotation for the ICGC cohort
#'
#' Built by `R/Load_PCAWG_clinical.R`; see the header of that file for what the
#' fields do and do not contain - in particular there is NO survival data, and
#' `project` doubles as a coarse receptor-status label that is confounded with
#' `grade`.
#'
#' @param donors Optional donor ids to select and order by, so the result lines
#'   up row-for-row with `colnames(fit$Thetas)` without a second `match()`.
#' @return A data frame, one row per donor.
load_clinical <- function(donors = NULL) {
  if (!file.exists(PATH_CLINICAL)) {
    stop("no clinical table at ", PATH_CLINICAL,
         "\n  Build it with:  Rscript R/Load_PCAWG_clinical.R", call. = FALSE)
  }
  out <- utils::read.csv(PATH_CLINICAL, stringsAsFactors = FALSE)
  if (is.null(donors)) return(out)

  gone <- setdiff(donors, out$donor)
  if (length(gone)) {
    stop(length(gone), " donor(s) have no clinical row: ",
         paste(utils::head(gone, 5), collapse = ", "), call. = FALSE)
  }
  out[match(donors, out$donor), , drop = FALSE]
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


#' Renumber a fit's signatures
#'
#' By default `SigN01` becomes the signature carrying the most mass, `SigN02` the
#' next, and so on, with every signature-indexed slot of the solution permuted to
#' agree. A fit numbers its signatures in whatever order the initialisation
#' happened to put them, which means the same process carries a different label
#' in every figure unless something imposes an order. This is that something, and
#' it is used by the de novo figures and the TensorSignatures comparison.
#'
#' Reorder the FIT rather than the matrices pulled out of it. Several things read
#' the labels back out of the fit itself and cannot see a reordered extract -
#' `plot_Mu()` rebuilds its panels from `df_assign(fit, data)`, and `df_assign()`
#' in turn reads `fit$Baseline`, NOT `fit$Thetas`, because the assignment
#' compares per-unit-exposure rates and \eqn{q_j(\beta_k)} would not cancel. A
#' hand-rolled permutation that moves `Signatures`, `Betas`, `Thetas` and `Mu`
#' but forgets `Baseline` still renders a perfectly plausible figure with every
#' mutation attributed to the wrong signature. This function moves all of them.
#'
#' Nothing changes numerically, but the rename is not safe to apply to half a
#' pipeline: a relabelled fit no longer agrees with a `df_assign()` table or a
#' figure built from the fit before it. Relabel once, straight after loading, and
#' pass the relabelled fit everywhere downstream. With the default order,
#' applying it twice is a no-op.
#'
#' @param fit A `SignaturePPF` fit.
#' @param order The new order, as positions into the fit's current signatures or
#'   as their names. `NULL` (default) sorts by decreasing `Mu`, which is what the
#'   function is named after. An explicit order is for figures that group
#'   signatures by what they are rather than by how big they are - putting SBS2
#'   next to SBS13, say - and is NOT idempotent, so apply it once.
#' @param prefix Label stem.
#' @return `fit`, with a `relabel` attribute recording old label -> new label.
relabel_by_mu <- function(fit, order = NULL, prefix = "SigN") {
  mu  <- as.numeric(fit$Mu)
  old <- colnames(fit$Signatures)

  if (is.null(order)) {
    ord <- base::order(mu, decreasing = TRUE)
  } else {
    ord <- if (is.character(order)) match(order, old) else as.integer(order)
    # A permutation that is silently wrong scrambles the fit without erroring, so
    # it is checked rather than trusted: every signature exactly once.
    if (anyNA(ord) || !setequal(ord, seq_along(old))) {
      stop("`order` must be a permutation of the fit's ", length(old),
           " signatures, given as positions or as names. Got: ",
           paste(utils::head(order, 12), collapse = ", "), call. = FALSE)
    }
  }
  new <- sprintf("%s%02d", prefix, seq_along(ord))

  perm_col <- function(m) { m <- m[, ord, drop = FALSE]; colnames(m) <- new; m }
  perm_row <- function(m) { m <- m[ord, , drop = FALSE]; rownames(m) <- new; m }

  for (s in intersect(c("Signatures", "SigPrior", "Betas"), names(fit)))
    fit[[s]] <- perm_col(fit[[s]])                     # I x K, I x K, p x K
  for (s in intersect(c("Thetas", "Baseline", "Q"), names(fit)))
    fit[[s]] <- perm_row(fit[[s]])                     # K x J
  for (s in intersect(c("Mu", "Sigma2"), names(fit)))
    fit[[s]] <- stats::setNames(as.numeric(fit[[s]])[ord], new)

  # The credible bounds are shaped like the block they bound.
  for (bound in c("lowCI", "highCI")) {
    ci <- fit[[bound]]
    if (is.null(ci)) next
    for (b in intersect(c("Signatures", "Betas"), names(ci))) ci[[b]] <- perm_col(ci[[b]])
    for (b in intersect("Thetas", names(ci))) ci[[b]] <- perm_row(ci[[b]])
    for (b in intersect(c("Mu", "Sigma2"), names(ci)))
      ci[[b]] <- stats::setNames(as.numeric(ci[[b]])[ord], new)
    fit[[bound]] <- ci
  }

  # The raw optimizer output carries no dimnames at all, so it cannot be caught
  # by a name lookup later: left alone it would silently disagree with the
  # relabelled solution above. It is only permutable while it still describes the
  # SAME signature set - after prune_signatures() the record deliberately keeps
  # every signature and the solution does not, and applying a K-long permutation
  # to a K_fitted-wide matrix would quietly scramble it.
  m <- fit$MAPsolution
  if (!is.null(m) && ncol(m$R) == length(ord)) {
    for (s in intersect(c("R", "Betas", "Mu", "Sigma2"), names(m)))
      m[[s]] <- m[[s]][, ord, drop = FALSE]            # Mu, Sigma2 are 1 x K
    for (s in intersect(c("Theta", "Phi", "Q"), names(m)))
      m[[s]] <- m[[s]][ord, , drop = FALSE]
    fit$MAPsolution <- m
  }

  # The chain DOES carry signature names, and `posterior_CI()` subsets it by the
  # solution's colnames - so leaving it under the old numbering would not error,
  # it would quietly return a different signature's draws. Relabel it to match.
  # Only the labels move: the draws stay in the order they were sampled, and any
  # signature `prune_signatures()` dropped is still there, numbered after the
  # retained block rather than renamed onto one of them.
  if (!is.null(fit$MCMCchain)) {
    full <- dimnames(fit$MCMCchain$MUchain)[[2]]
    lab <- stats::setNames(rep(NA_character_, length(full)), full)
    lab[old[ord]] <- new
    rest <- full[is.na(lab)]
    if (length(rest)) {
      lab[rest] <- sprintf("%s%02d", prefix, length(ord) + seq_along(rest))
    }
    ch <- fit$MCMCchain
    for (s in c("SIGSchain", "BETASchain")) {                 # draw x . x sig
      dimnames(ch[[s]])[[3]] <- unname(lab[dimnames(ch[[s]])[[3]]])
    }
    for (s in c("THETAchain")) {                              # draw x sig x sample
      dimnames(ch[[s]])[[2]] <- unname(lab[dimnames(ch[[s]])[[2]]])
    }
    for (s in c("MUchain", "SIGMA2chain", "SHRINKchain", "ADAPTchain")) {
      if (is.null(ch[[s]])) next
      dimnames(ch[[s]])[[2]] <- unname(lab[dimnames(ch[[s]])[[2]]])
    }
    fit$MCMCchain <- ch
    if (length(fit$pruned)) fit$pruned <- unname(lab[fit$pruned])
  }

  attr(fit, "relabel") <- data.frame(from = old[ord], to = new, mu = mu[ord],
                                     row.names = NULL, stringsAsFactors = FALSE)
  fit
}



## Mutation assignment lives in the package now: assign_mutations(), df_assign()
## and the phi lookup they share are all exported by SignaturePPF. Keeping copies
## here masked them whenever load_functions() ran after library(SignaturePPF),
## which is every script - and a masked copy is worse than a missing one, because
## it silently wins.


## ---------------------------------------------------------------------------
## MCMC bookkeeping shared by the application figure scripts and the simulation
## studies. These lived in Simulation_functions_misspec.R, which the application
## scripts never source, so Reproduce_figures_Application_denovo.R and _refit.R
## both failed on them. They are not misspecification code, so they belong here,
## where load_functions() picks them up for every script.
## ---------------------------------------------------------------------------

#' Which stored draws survive the burn-in
#'
#' The chain is indexed by stored draw and `burnin` is in iterations, so the cut is
#' at `index * thin > burnin`. Dropping the first `burnin` ROWS would discard the
#' wrong draws whenever `thin > 1`.
kept_draw_index <- function(fit) {
  n <- dim(fit$MCMCchain$MUchain)[1]
  which(seq_len(n) * fit$controls$thin > fit$controls$burnin)
}

#' Effective sample size of a chain block, post burn-in
get_PosteriorEffectiveSize <- function(chain, keep) {
  # The log-density chains are only recorded every `logpost_every` iterations,
  # so the slots in between are NA, and coda::effectiveSize calls na.fail on
  # them. Dropping them makes the ESS that of the recorded (thinned) series,
  # which is the only series that exists. The parameter chains are dense, so
  # this is a no-op for them.
  ess1 <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) NA_real_ else unname(coda::effectiveSize(x))
  }
  if (is.null(dim(chain))) {
    ess1(chain[keep])
  } else if (length(dim(chain)) == 2) {
    apply(chain[keep, , drop = FALSE], 2, ess1)
  } else {
    apply(chain[keep, , , drop = FALSE], c(2, 3), ess1)
  }
}

#' Attribution probabilities of every mutation under a fitted model
#'
#' `Theta` here is the BASELINE. `Betas = NULL` gives the position-independent
#' NMF case, where the probability reduces to R[i,k] * theta[k,j].
Compute_mutation_Probs <- function(gr_Mutations, R, Theta, Betas = NULL) {
  ch <- as.character(gr_Mutations$channel)
  sm <- as.character(gr_Mutations$sample)
  Probs <- R[ch, , drop = FALSE] * t(Theta[, sm, drop = FALSE])
  if (!is.null(Betas)) {
    X <- as.matrix(GenomicRanges::mcols(gr_Mutations)[, rownames(Betas), drop = FALSE])
    Probs <- Probs * exp(X %*% Betas)
  }
  Probs <- Probs / rowSums(Probs)
  colnames(Probs) <- colnames(R)
  rownames(Probs) <- ch
  Probs
}
