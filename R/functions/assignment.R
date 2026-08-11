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
#' @param data The cohort it was fitted on.
#' @return A character vector, one signature name per mutation.
assign_mutations <- function(fit, data) {
  gr <- data$gr_Mutations
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
