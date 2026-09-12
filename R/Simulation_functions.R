## Produces: no figure. Helper file, sourced by the analysis scripts.

## Generative core of the simulation study.
##
## Carried over from the predecessor project essentially unchanged - these
## functions define the data-generating process the paper reports, so they are
## deliberately NOT rewritten. Three edits only, each marked at its site:
##
##   * `SigPoisProcess::` -> `SignaturePPF::` for the COSMIC catalogue;
##   * `generate_Parameters(K_new = 0)` no longer references an undefined
##     `Rmat_random`;
##   * `cosine()` is qualified as `sigminer::cosine()`, so the file works when
##     sourced rather than run after `library(sigminer)`.
##
## The misspecification variants live in Simulation_functions_misspec.R.

#' Read an RDS if it exists, otherwise NULL
open_rds_file <- function(file) {
  if (file.exists(file)) readRDS(file) else NULL
}

#' Create a directory if it is not already there
create_directory <- function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}


################################################################################
# Part 1 - Simulating the data
################################################################################

#' Genomic covariate track: an AR(1) Gaussian process along a single chromosome
#'
#' `rho` close to 1 makes neighbouring tiles strongly correlated, which is what
#' real epigenomic tracks look like and what makes the covariates hard to
#' separate. Columns are standardised, matching the real preprocessing.
generate_SignalTrack <- function(length_genome = 2e6,
                                 tilewidth = 100,
                                 rho = 0.99, p = 10,
                                 Sigma = diag(p)) {
  grX <- GenomicRanges::tileGenome(c("chrsim" = length_genome),
                                   tilewidth = tilewidth,
                                   cut.last.tile.in.chrom = TRUE)
  X <- matrix(NA, nrow = length(grX), ncol = p)
  D <- diag(rep(rho, p))
  X[1, ] <- mvtnorm::rmvnorm(1, sigma = Sigma)
  for (i in 2:nrow(X)) {
    eps <- c(mvtnorm::rmvnorm(1, sigma = Sigma))
    X[i, ] <- c(D %*% X[i - 1, ]) + c(eps)
  }
  colnames(X) <- paste0("Encode", sprintf("%02d", 1:ncol(X)))
  grX$bin_weight <- tilewidth
  GenomicRanges::mcols(grX) <- cbind(GenomicRanges::mcols(grX), scale(X))
  grX
}


#' Per-patient copy number: piecewise-constant segments of gains
generate_CopyTrack <- function(J = 50,
                               mu_copy = 1,
                               size_copy = 1,
                               length_genome = 2e6,
                               tilewidth = 100,
                               mean_segments = 20) {
  grX <- GenomicRanges::tileGenome(c("chrsim" = length_genome),
                                   tilewidth = tilewidth,
                                   cut.last.tile.in.chrom = TRUE)
  n_tiles <- length(grX)
  cnmat <- matrix(NA, nrow = n_tiles, ncol = J)
  for (j in seq_len(J)) {
    n_seg <- rpois(1, lambda = mean_segments)
    segs <- sort(sample(1:n_tiles, n_seg, replace = FALSE))
    segs <- c(1, segs, n_tiles + 1)
    profile <- numeric(n_tiles)
    for (k in seq_len(length(segs) - 1)) {
      s <- segs[k]
      e <- segs[k + 1] - 1
      if (e < s) next
      cn <- 2 + rnbinom(1, mu = mu_copy, size = size_copy)
      profile[s:e] <- cn
    }
    cnmat[, j] <- profile
  }
  colnames(cnmat) <- paste0("Patient_", sprintf("%02d", 1:J))
  GenomicRanges::mcols(grX) <- as.data.frame(cnmat)
  grX
}


#' True parameters: signatures R, activities Theta and coefficients Betas
#'
#' `p_all - p_to_use` covariates are given a coefficient of exactly zero. They
#' are the redundant covariates the false-positive rate is measured on.
generate_Parameters <- function(cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                K_new = 4,
                                J = 50,
                                a = 1,
                                p_all = 10,
                                p_to_use = 5,
                                theta = 200,
                                prob_zero = 0,
                                sd_beta = 1 / 2,
                                alpha = 0.1) {
  p <- p_to_use
  Rmat_cos <- SignaturePPF::COSMIC_v3.4_SBS96_GRCh37[, cosmic_sigs]

  # FIXED: the original unconditionally did cbind(Rmat_cos, Rmat_random), which
  # errors with "object 'Rmat_random' not found" whenever K_new = 0 - the setting
  # the misspecification study is run in. See the note in the study script.
  if (K_new > 0) {
    Rmat_random <- t(LaplacesDemon::rdirichlet(n = K_new, alpha = rep(alpha, 96)))
    colnames(Rmat_random) <- paste0("SBSnew", 1:K_new)
    Rmat <- cbind(Rmat_cos, Rmat_random)
  } else {
    Rmat <- Rmat_cos
  }

  K_true <- ncol(Rmat)
  exposures <- rgamma(K_true, theta, 1)
  Theta <- matrix(rgamma(K_true * J, 0.5, 0.5), ncol = J, nrow = K_true)
  Theta <- apply(Theta, 2, function(x) x * exposures)
  rownames(Theta) <- colnames(Rmat)

  Betas <- matrix(rnorm(K_true * p, mean = 0, sd = sd_beta), nrow = p, ncol = K_true)
  Null <- 1 - matrix(rbinom(n = length(Betas), size = 1, prob = prob_zero),
                     nrow = nrow(Betas), ncol = ncol(Betas))
  if (any(colSums(Null) == 0)) {
    id <- sample(1:p, size = 1)
    Null[id, colSums(Null) == 0] <- 1
  }
  Betas <- Betas * Null
  if (p_all > p_to_use) {
    Betas <- rbind(Betas, matrix(0, nrow = p_all - p_to_use, ncol = K_true))
  }

  list(R = Rmat, Theta = Theta, Betas = Betas)
}


#' Expected counts per bin and patient under a set of parameters
#'
#' Replaces `SigPoisProcess:::Reconstruct_Lambda`, which SignaturePPF does not
#' carry. The C++ version computed
#' \eqn{\Lambda_{tj} = \sum_i \sum_k c_j(t)\,e^{\beta_k'x(t)}\phi_{kj} r_{ik}};
#' the sum over channels collapses because the columns of \eqn{R} sum to one, so
#' this is one matrix product and agrees with the old code to machine precision.
#'
#' NOTE the third argument is the BASELINE \eqn{\phi}, the per-unit-exposure
#' activity - not the total. Copy number here already carries the tile width and
#' the factor of one half, as the callers set it up.
reconstruct_lambda_raw <- function(SignalTrack, CopyTrack, Phi, Betas) {
  E <- exp(pmin(pmax(SignalTrack[, rownames(Betas), drop = FALSE] %*% Betas,
                     -20), 20))
  CopyTrack * (E %*% Phi)
}


################################################################################
# Part 2 - Scoring the output
################################################################################

#' Sensitivity / precision / F1 of a recovered signature set
Compute_sensitivity_precision <- function(R_hat, R_true, cos_cutoff = 0.9) {
  sig_sens <- sapply(1:ncol(R_true), function(i)
    max(sapply(1:ncol(R_hat), function(j) sigminer::cosine(R_true[, i], R_hat[, j]))))
  sig_prec <- sapply(1:ncol(R_hat), function(i)
    max(sapply(1:ncol(R_true), function(j) sigminer::cosine(R_true[, j], R_hat[, i]))))
  sens <- mean(sig_sens > cos_cutoff)
  prec <- mean(sig_prec > cos_cutoff)
  f1 <- if (sens + prec > 0) 2 * sens * prec / (sens + prec) else 0
  c("Sensitivity" = sens, "Precision" = prec, "F1" = f1)
}


#' One-to-one matching of estimated to true signatures
#'
#' The smaller matrix is padded with a CONSTANT (100) rather than with zeros: a
#' zero column has no direction, so its cosine against anything is NaN and the
#' assignment problem becomes unsolvable.
match_MutSign <- function(R_true, R_hat) {
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  k_tot <- max(c(k_hat, k_true))
  I <- nrow(R_true)
  mat0 <- matrix(100, nrow = I, ncol = abs(k_hat - k_true))
  if (k_hat > k_true) {
    colnames(mat0) <- paste0("new_extra", 1:ncol(mat0))
    R_true <- cbind(R_true, mat0)
  } else if (k_hat < k_true) {
    R_hat <- cbind(R_hat, mat0)
  }

  CosMat <- matrix(1, k_tot, k_tot)
  for (i in 1:k_tot) {
    for (j in 1:k_tot) {
      CosMat[i, j] <- 1 - sigminer::cosine(R_true[, i], R_hat[, j])
    }
  }
  match <- RcppHungarian::HungarianSolver(CosMat)$pairs[, 2]
  R_hat_matched <- R_hat[, match]

  R_hat_matched[R_hat_matched == 100] <- 0
  R_true[R_true == 100] <- 0

  list("R_hat" = R_hat_matched, "R_true" = R_true, "match" = match)
}


#' RMSE of the activities, after matching
compute_RMSE_Theta <- function(Theta_true, Theta_hat, match) {
  k_true <- nrow(Theta_true)
  k_hat <- nrow(Theta_hat)
  J <- ncol(Theta_true)
  mat0 <- matrix(0, nrow = abs(k_hat - k_true), ncol = J)
  if (k_hat > k_true) {
    rownames(mat0) <- paste0("new_extra", 1:nrow(mat0))
    Theta_true <- rbind(Theta_true, mat0)
  } else if (k_hat < k_true) {
    Theta_hat <- rbind(Theta_hat, mat0)
  }
  Theta_hat <- Theta_hat[match, ]
  sqrt(mean((Theta_true - Theta_hat)^2))
}


#' RMSE of the covariate coefficients, after matching
compute_RMSE_Betas <- function(Beta_true, Beta_hat, match) {
  p_true <- nrow(Beta_true)
  k_true <- ncol(Beta_true)
  p_hat <- nrow(Beta_hat)
  k_hat <- ncol(Beta_hat)
  p_max <- max(p_true, p_hat)
  k_max <- max(k_true, k_hat)
  pad_matrix <- function(mat, nrows, ncols) {
    out <- matrix(0, nrows, ncols)
    out[1:nrow(mat), 1:ncol(mat)] <- mat
    out
  }
  Beta_true_pad <- pad_matrix(Beta_true, p_max, k_max)
  Beta_hat_pad <- pad_matrix(Beta_hat, p_max, k_max)
  Beta_hat_matched <- Beta_hat_pad[, match, drop = FALSE]
  Beta_true_matched <- Beta_true_pad[, seq_along(match), drop = FALSE]
  stopifnot(all(dim(Beta_true_matched) == dim(Beta_hat_matched)))
  sqrt(mean((Beta_true_matched - Beta_hat_matched)^2))
}


#' Reshape a matrix onto the dimensions and dimnames of a reference
pad_to_match <- function(x, reference, fill = 0, strict = FALSE) {
  x <- as.matrix(x)
  ref <- as.matrix(reference)

  if (strict) {
    extra_r <- setdiff(rownames(x), rownames(ref))
    extra_c <- setdiff(colnames(x), colnames(ref))
    if (length(extra_r)) stop("rows in x not in reference: ", paste(extra_r, collapse = ", "))
    if (length(extra_c)) stop("cols in x not in reference: ", paste(extra_c, collapse = ", "))
  }

  out <- matrix(fill, nrow = nrow(ref), ncol = ncol(ref), dimnames = dimnames(ref))
  common_r <- intersect(rownames(ref), rownames(x))
  common_c <- intersect(colnames(ref), colnames(x))
  out[common_r, common_c] <- x[common_r, common_c]
  out
}


#' Add zero columns so a matrix carries a required column set, in order
add_missing_cols <- function(x, all_names, fill = 0) {
  missing <- setdiff(all_names, colnames(x))
  if (length(missing)) {
    zeros <- matrix(fill, nrow = nrow(x), ncol = length(missing),
                    dimnames = list(rownames(x), missing))
    x <- cbind(x, zeros)
  }
  x[, all_names, drop = FALSE]
}
