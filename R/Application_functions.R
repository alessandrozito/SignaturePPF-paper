## Shared machinery for the three ICGC applications.
##
## All of them follow the same shape: find a mode, then sample from it. What
## differs is how many starting points the mode is searched from, and whether the
## signatures are estimated or held fixed.

#' Log posterior at the end of a MAP run
#'
#' The optimiser's objective trace, whose last entry is the value it stopped at.
#' Evaluated every ten iterations rather than every one, so the trace is shorter
#' than the iteration count - that is a diagnostic subsample, not a thinned chain.
map_logposterior <- function(fit) {
  tr <- fit$MAPsolution$trace
  if (is.null(tr) || !length(tr)) return(NA_real_)
  utils::tail(as.numeric(tr), 1)
}


#' MAP from several random starts, keeping the best
#'
#' The de novo posterior is multimodal - signatures can swap, split and merge -
#' so a single optimisation reports whichever basin its starting point fell into.
#' Running a few and keeping the highest log posterior is the cheap defence, and
#' it also makes the spread across starts visible, which is the honest thing to
#' report: if the three modes are far apart, that is worth knowing.
#'
#' Every start is cached separately, so an interrupted search resumes and a
#' fourth start can be added later without redoing the first three.
#'
#' @param data A cohort object.
#' @param out_dir Where the per-start fits and the summary table are written.
#' @param n_starts Number of random starting points.
#' @param seed Base seed. Start `i` uses `seed + i`, so the set is reproducible
#'   and extending it does not change the earlier starts.
#' @param ... Passed to [SignaturePPF::SignaturePPF()] - `sigs`, `K`, `prior`,
#'   `controls` and so on.
#' @return The fit with the highest log posterior, with `$start` recording which
#'   one it was.
fit_map_restarts <- function(data, out_dir, n_starts = 3, seed = SEED,
                             verbose = TRUE, ...) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fits <- lapply(seq_len(n_starts), function(i) {
    file <- file.path(out_dir, sprintf("MAPSolution_start%02d.rds.gzip", i))
    if (file.exists(file)) {
      message("using existing MAP start ", i, ": ", basename(file))
      return(readRDS(file))
    }
    message("\n=== MAP start ", i, " of ", n_starts, " ===")
    fit <- SignaturePPF::SignaturePPF(data, method = "map", seed = seed + i,
                                      verbose = verbose, ...)
    saveRDS(fit, file, compress = "gzip")
    fit
  })

  lp <- vapply(fits, map_logposterior, numeric(1))
  summary_tbl <- data.frame(
    start = seq_len(n_starts),
    logposterior = lp,
    iterations = vapply(fits, function(f) as.integer(f$MAPsolution$it), integer(1)),
    n_active = vapply(fits, function(f) sum(f$Mu > 10 * f$prior$epsilon), integer(1)),
    minutes = vapply(fits, function(f) as.numeric(f$runtime, units = "mins"), numeric(1)))
  summary_tbl$best <- seq_len(n_starts) == which.max(lp)
  utils::write.csv(summary_tbl, file.path(out_dir, "map_starts.csv"),
                   row.names = FALSE)
  print(summary_tbl)

  best <- which.max(lp)
  message("\nkeeping start ", best, " (log posterior ",
          format(lp[best], big.mark = ",", nsmall = 2), ")")
  out <- fits[[best]]
  out$start <- best
  out
}


#' Initial values for an MCMC chain, taken from a MAP fit
#'
#' `R_start` is omitted for a refit: the signatures come from `sigs` there, and
#' passing both is rejected by the package rather than quietly ignored.
#'
#' @param fit A MAP fit.
#' @param sigs_fixed Whether the chain that follows holds the signatures fixed.
init_from_map <- function(fit, sigs_fixed = FALSE) {
  args <- list(Theta_start = fit$Thetas,     # activity scale, as the package expects
               Betas_start = fit$Betas,
               Mu_start = fit$Mu,
               Sigma2_start = fit$Sigma2)
  if (!sigs_fixed) args$R_start <- fit$Signatures
  do.call(SignaturePPF::SignaturePPF_init, args)
}


#' Run the chain from a MAP fit, checkpointing as it goes
#'
#' The predecessor did this by hand: run 100 sweeps, `abind` them onto the
#' chain so far, re-save the whole growing object, repeat. That rewrote the
#' entire chain to disk every hundred sweeps and could not restore the RNG, so a
#' resumed run was a different chain that happened to start where the old one
#' stopped. `SignaturePPF_checkpoint()` writes each block once and restores
#' `.Random.seed` with it, so an interrupted run resumes bit-exactly.
#'
#' @param data A cohort object.
#' @param map A MAP fit to start from.
#' @param out_dir Where the chain and its checkpoint directory live.
#' @param sigs Reference signatures, for a refit. `NULL` for de novo.
#' @param every Sweeps between checkpoints.
#' @param ... Passed to [SignaturePPF::SignaturePPF()].
#' @return The fitted chain.
run_mcmc_from_map <- function(data, map, out_dir, sigs = NULL,
                              sigs_fixed = !is.null(sigs), every = 100,
                              seed = SEED, verbose = TRUE, ...) {
  out_file <- file.path(out_dir, "MCMCSolution.rds.gzip")
  if (file.exists(out_file)) {
    message("using existing chain: ", basename(out_file))
    return(readRDS(out_file))
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- SignaturePPF::SignaturePPF(
    data,
    sigs = sigs,
    sigs_fixed = sigs_fixed,
    method = "mcmc",
    init = init_from_map(map, sigs_fixed = sigs_fixed),
    # The mode has already been found and paid for; re-running the optimiser
    # inside the MCMC call would just find it again.
    init_mcmc_from_map = FALSE,
    # Signatures the compressive prior has parked near epsilon are KEPT. Pruning
    # them would fix the model's dimension at the mode's answer, and whether a
    # marginal signature survives is one of the things the posterior is being
    # asked about.
    prune_after_map = FALSE,
    checkpoint = SignaturePPF::SignaturePPF_checkpoint(
      dir = file.path(out_dir, "checkpoint"), every = every, resume = TRUE),
    seed = seed,
    verbose = verbose,
    ...)

  saveRDS(fit, out_file, compress = "gzip")
  fit
}


#' Posterior summaries of a chain, in the shape the figures expect
#'
#' @param fit An MCMC fit.
#' @param level Credible level.
posterior_summaries <- function(fit, level = 0.95) {
  blocks <- c("Signatures", "Thetas", "Betas", "Mu", "Sigma2")
  out <- lapply(blocks, function(b) SignaturePPF::posterior_CI(fit, b, level))
  names(out) <- blocks
  out$logPost <- fit$MCMCchain$logPostchain
  out$logLik <- fit$MCMCchain$logLikchain
  out$logPrior <- fit$MCMCchain$logPriorchain
  out$n_kept <- fit$n_kept
  out$level <- level
  out
}


#' Trace of the log posterior, with the burn-in marked
#'
#' The chain is indexed by STORED DRAW, so the sweep a draw corresponds to is its
#' index times `thin` - and the burn-in is specified in sweeps. Plotting against
#' the raw index would put the burn-in line in the wrong place by that factor.
#' `logpost_every` does not enter: it leaves NA in the entries where the
#' objective was not evaluated, and those are simply dropped.
plot_logposterior_trace <- function(fit) {
  lp <- fit$MCMCchain$logPostchain
  keep <- is.finite(lp)
  df <- data.frame(sweep = which(keep) * fit$controls$thin, logpost = lp[keep])

  ggplot2::ggplot(df, ggplot2::aes(.data$sweep, .data$logpost)) +
    ggplot2::geom_vline(xintercept = fit$controls$burnin, linetype = 2,
                        colour = "#CD2626") +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::labs(x = "Sweep", y = "Log posterior") +
    ggplot2::theme_bw()
}
