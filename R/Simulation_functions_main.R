## The main simulation study (Section 4): data generation, fitting and scoring.
##
## Companion to Simulation_functions.R, which holds the generative core shared
## with the misspecification study. As there, the generative half is carried over
## unchanged and the fitting/scoring half is ported to SignaturePPF - where
## `$Thetas` is the TOTAL activity and the per-unit-exposure baseline is
## `$Baseline`.

################################################################################
# Part 1 - Generating the data
################################################################################

#' Mutations of one channel for one patient
#'
#' Positions are drawn by inverting the cumulative intensity, so the points land
#' with the right spatial distribution rather than being binned after the fact.
sample_mutations <- function(i, j, R, Theta, RangesX, ExpBetas, CopyTrack,
                             tilewidth = 100) {
  CumSums <- apply(ExpBetas, 2, function(x) cumsum(0.5 * tilewidth * CopyTrack[, j] * x))
  Lambda_Tmax <- c(crossprod(R[i, ], Theta[, j] * c(utils::tail(CumSums, 1))))
  Lambdas <- c(crossprod(R[i, ], Theta[, j] * t(CumSums)))
  N_Tmax <- rpois(1, Lambda_Tmax)
  if (N_Tmax > 0) {
    u <- sort(runif(n = N_Tmax))
    findInterval(u, Lambdas / Lambda_Tmax) * tilewidth +
      sample(1:(tilewidth - 1), size = length(u), replace = TRUE)
  } else {
    NA
  }
}


#' Every patient's mutations
#'
#' As in the misspecification study, `ncores = 1` registers the sequential
#' backend: the driver parallelises over replicates, so a nested fork per patient
#' would oversubscribe the machine.
sample_dataset_parallelCN <- function(R, Theta, Betas, RangesX, ExpBetas,
                                      CopyTrack, tilewidth, ncores = 1,
                                      verbose = FALSE) {
  I <- nrow(R)
  J <- ncol(Theta)

  if (ncores > 1) doParallel::registerDoParallel(ncores) else foreach::registerDoSEQ()
  j <- NULL
  Mutations <- foreach::foreach(j = 1:J, .combine = "rbind") %dopar% {
    if (verbose) message("Simulating Patient_", sprintf("%02d", j))
    mut_temp <- data.frame()
    for (i in 1:I) {
      ch <- rownames(R)[i]
      ts <- sample_mutations(i, j, R, Theta, RangesX, ExpBetas, CopyTrack, tilewidth)
      if (!is.na(ts[1])) {
        mut_temp <- rbind(mut_temp,
                          data.frame(sample = paste0("Patient_", sprintf("%02d", j)),
                                     pos = ts, channel = ch))
      }
    }
    mut_temp
  }
  Mutations$sample <- factor(Mutations$sample, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)
  Mutations
}


#' One dataset for the main simulation
#'
#' `p_to_use` of the `p_all` covariates actually drive the intensity; the rest
#' have coefficients of exactly zero and are the redundant covariates the study
#' asks whether the model can decline.
generate_MutationData <- function(J = 50,
                                  cosmic_sigs = c("SBS1", "SBS2", "SBS13", "SBS3"),
                                  K_new = 4,
                                  theta = 100,
                                  a = 1,
                                  mu_copy = 1,
                                  size_copy = 10,
                                  sd_beta = 1 / 2,
                                  length_genome = 2e6,
                                  tilewidth = 100,
                                  rho = 0.99,
                                  p_all = 10,
                                  p_to_use = 5,
                                  corr = "indep",
                                  ncores = 1) {

  if (corr == "indep") {
    Sigma <- diag(p_all)
  } else if (corr == "onion") {
    corrmat <- clusterGeneration::genPositiveDefMat(p_all, covMethod = "onion")$Sigma
    Sigma <- cov2cor(corrmat)
  }

  grX <- generate_SignalTrack(length_genome = length_genome,
                              tilewidth = tilewidth, rho = rho, p = p_all,
                              Sigma = Sigma)
  gr_CopyTrack <- generate_CopyTrack(J = J, mu_copy = mu_copy,
                                     size_copy = size_copy,
                                     length_genome = length_genome,
                                     tilewidth = tilewidth)
  pars_true <- generate_Parameters(cosmic_sigs = cosmic_sigs, K_new = K_new,
                                   J = J, a = a, p_all = p_all,
                                   p_to_use = p_to_use, theta = theta,
                                   sd_beta = sd_beta)

  Xcovs <- as.matrix(GenomicRanges::mcols(grX)[, -1])   # drop bin_weight
  bin_weight <- grX$bin_weight
  ExpBetas <- exp(Xcovs %*% pars_true$Betas)
  CopyTrack <- as.matrix(GenomicRanges::mcols(gr_CopyTrack))
  colnames(pars_true$Theta) <- colnames(CopyTrack)
  Theta_scaled <- pars_true$Theta / sum(bin_weight)

  Mutations <- sample_dataset_parallelCN(R = pars_true$R,
                                         Betas = pars_true$Betas,
                                         Theta = Theta_scaled,
                                         RangesX = IRanges::ranges(grX),
                                         ExpBetas = ExpBetas,
                                         CopyTrack = CopyTrack,
                                         tilewidth = tilewidth,
                                         ncores = ncores)

  gr_Mutations <- GenomicRanges::GRanges(
    seqnames = "chrsim",
    IRanges::IRanges(start = Mutations$pos, end = Mutations$pos),
    sample = Mutations$sample, channel = Mutations$channel)
  overlaps <- GenomicRanges::findOverlaps(gr_Mutations, grX)
  GenomicRanges::mcols(gr_Mutations) <- cbind(
    GenomicRanges::mcols(gr_Mutations),
    GenomicRanges::mcols(grX[S4Vectors::subjectHits(overlaps)]))
  gr_Mutations$bin_weight <- NULL

  Betas <- pars_true$Betas
  rownames(Betas) <- colnames(Xcovs)
  colnames(Betas) <- colnames(pars_true$R)

  list("gr_Mutations" = gr_Mutations,
       "covariate_used" = colnames(Xcovs)[1:p_to_use],
       "R" = pars_true$R,
       "Betas" = Betas,
       "Theta" = Theta_scaled,
       "gr_SignalTrack" = grX,
       "gr_CopyTrack" = gr_CopyTrack,
       "simulation_parameters" = list(J = J, cosmic_sigs = cosmic_sigs,
                                      K_new = K_new, theta = theta, a = a,
                                      mu_copy = mu_copy, size_copy = size_copy,
                                      length_genome = length_genome,
                                      tilewidth = tilewidth, rho = rho,
                                      p_all = p_all, p_to_use = p_to_use,
                                      corr = corr))
}


#' Coarsen the covariates and copy number onto wider bins
#'
#' Models M5 and M6 are fitted to covariates averaged over groups of consecutive
#' tiles, which is what a real analysis faces: the tracks are measured at one
#' resolution and the mutations at another.
#'
#' Two versions of each track are returned. `SignalTrack` / `CopyTrack` are the
#' COARSE grid, one row per wide bin, and are what the model is fitted on.
#' `SignalTrackRed` / `CopyTrackRed` are the same values broadcast back onto the
#' ORIGINAL fine grid, and are what the scoring uses - the truth lives on the
#' fine grid, so a coarse-grid intensity has to be expanded before it can be
#' compared with it.
aggregate_SignalTrack_CopyTrack <- function(data, lenght_bins = 200) {
  df <- as.data.frame(data$gr_SignalTrack)
  df_copy <- as.data.frame(data$gr_CopyTrack)

  df$bin_start <- df_copy$bin_start <-
    floor((df$start - 1) / lenght_bins) * lenght_bins + 1
  df$bin_end <- df_copy$bin_end <- df$bin_start + lenght_bins - 1

  df_agg_sum <- df |>
    dplyr::group_by(.data$seqnames, .data$bin_start, .data$bin_end) |>
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric), mean),
                     .groups = "drop")
  gr_agg <- GenomicRanges::GRanges(
    seqnames = df_agg_sum$seqnames,
    ranges = IRanges::IRanges(start = df_agg_sum$bin_start,
                              end = df_agg_sum$bin_end))
  GenomicRanges::mcols(gr_agg) <- df_agg_sum[, -(1:6)]
  gr_agg$bin_weight <- lenght_bins

  df_copy_sum <- df_copy |>
    dplyr::group_by(.data$seqnames, .data$bin_start, .data$bin_end) |>
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric), mean),
                     .groups = "drop")
  gr_copy <- GenomicRanges::GRanges(
    seqnames = df_copy_sum$seqnames,
    ranges = IRanges::IRanges(start = df_copy_sum$bin_start,
                              end = df_copy_sum$bin_end))
  GenomicRanges::mcols(gr_copy) <- df_copy_sum[, -(1:6)]

  # `mutate`, not `summarize`: the coarse value written back onto every fine bin.
  df_agg_tot <- df |>
    dplyr::group_by(.data$seqnames, .data$bin_start, .data$bin_end) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), mean))
  gr_agg_tot <- GenomicRanges::GRanges(
    seqnames = df_agg_tot$seqnames,
    ranges = IRanges::IRanges(start = df_agg_tot$bin_start,
                              end = df_agg_tot$bin_end))
  GenomicRanges::mcols(gr_agg_tot) <- df_agg_tot[, -(1:6)] |>
    dplyr::ungroup() |> dplyr::select(-"bin_start", -"bin_end")

  df_copy_tot <- df_copy |>
    dplyr::group_by(.data$seqnames, .data$bin_start, .data$bin_end) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), mean))
  gr_copy_tot <- GenomicRanges::GRanges(
    seqnames = df_copy_tot$seqnames,
    ranges = IRanges::IRanges(start = df_copy_tot$bin_start,
                              end = df_copy_tot$bin_end))
  GenomicRanges::mcols(gr_copy_tot) <- df_copy_tot[, -(1:5)] |>
    dplyr::ungroup() |> dplyr::select(-"bin_start", -"bin_end")

  gr_Mutations_agg <- data$gr_Mutations
  overlaps <- GenomicRanges::findOverlaps(gr_Mutations_agg, gr_agg)
  GenomicRanges::mcols(gr_Mutations_agg) <- cbind(
    GenomicRanges::mcols(gr_Mutations_agg)[, c(1, 2)],
    GenomicRanges::mcols(gr_agg[S4Vectors::subjectHits(overlaps)]))
  gr_Mutations_agg$bin_weight <- NULL

  list(gr_Mutations = gr_Mutations_agg,
       SignalTrack = as.matrix(GenomicRanges::mcols(gr_agg))[, -1],
       SignalTrackRed = as.matrix(GenomicRanges::mcols(gr_agg_tot)),
       CopyTrack = as.matrix(GenomicRanges::mcols(gr_copy)) / 2 * lenght_bins,
       CopyTrackRed = as.matrix(GenomicRanges::mcols(gr_copy_tot)) / 2 *
         data$simulation_parameters$tilewidth)
}


################################################################################
# Part 2 - Scoring
################################################################################

#' Observed counts per bin and sample
get_counts_from_data <- function(data) {
  over <- GenomicRanges::findOverlaps(data$gr_Mutations, data$gr_SignalTrack)
  region <- S4Vectors::subjectHits(over)
  tab <- table(factor(region, levels = seq_along(data$gr_SignalTrack)),
               factor(as.character(data$gr_Mutations$sample),
                      levels = levels(data$gr_Mutations$sample)))
  matrix(as.numeric(tab), nrow = length(data$gr_SignalTrack),
         dimnames = list(NULL, colnames(tab)))
}


#' Which stored draws survive the burn-in
#'
#' The chain is indexed by stored draw and `burnin` is in sweeps, so the cut is
#' at `index * thin > burnin`. Dropping the first `burnin` ROWS would discard the
#' wrong draws whenever `thin > 1`.
kept_draw_index <- function(fit) {
  n <- dim(fit$MCMCchain$MUchain)[1]
  which(seq_len(n) * fit$controls$thin > fit$controls$burnin)
}


#' Effective sample size of a chain block, post burn-in
get_PosteriorEffectiveSize <- function(chain, keep) {
  if (is.null(dim(chain))) {
    coda::effectiveSize(chain[keep])
  } else if (length(dim(chain)) == 2) {
    apply(chain[keep, , drop = FALSE], 2, coda::effectiveSize)
  } else {
    apply(chain[keep, , , drop = FALSE], c(2, 3), coda::effectiveSize)
  }
}


#' A common estimate structure for any of the fitted models
#'
#' `model = "CompNMF"` and `"SignatureAnalyzer"` are the covariate-free
#' competitors: their exposures are spread over the genome by copy number alone,
#' which is the finest prediction they can express.
getModelEstimates <- function(res, model = "SignaturePPF", SignalTrack, CopyTrack) {
  flat <- matrix(rep(1 / 96, 96))

  if (model == "CompNMF") {
    W <- res$Signatures; H <- res$Theta
    keep <- (res$Mu > 0) & c(sigminer::cosine(W, flat) < 0.975)
  } else if (model == "SignatureAnalyzer") {
    W <- res$Signature.norm; H <- res$Exposure
    keep <- c(sigminer::cosine(W, flat) < 0.975)
  } else {
    cutoff <- 5 * res$prior$a * res$prior$epsilon
    keep <- (res$Mu > 5 * cutoff) & c(sigminer::cosine(res$Signatures, flat) < 0.975)
  }
  if (!any(keep)) keep <- rep(TRUE, length(keep))

  if (model %in% c("CompNMF", "SignatureAnalyzer")) {
    R_hat <- W[, keep, drop = FALSE]
    Theta_total <- H[keep, , drop = FALSE]
    theta_baseline <- Theta_total /
      t(colSums(CopyTrack))[rep(1, nrow(Theta_total)), , drop = FALSE]
    Beta_hat <- matrix(0, nrow = ncol(SignalTrack), ncol = nrow(theta_baseline),
                       dimnames = list(colnames(SignalTrack), rownames(theta_baseline)))
  } else {
    R_hat <- res$Signatures[, keep, drop = FALSE]
    Beta_hat <- res$Betas[, keep, drop = FALSE]
    # $Baseline, NOT $Thetas: under the activity prior the latter is the total.
    theta_baseline <- res$Baseline[keep, , drop = FALSE]
    Theta_total <- theta_baseline *
      crossprod(exp(SignalTrack[, rownames(Beta_hat), drop = FALSE] %*% Beta_hat),
                CopyTrack)
  }
  Lambda_hat <- reconstruct_lambda_raw(SignalTrack = SignalTrack,
                                       CopyTrack = CopyTrack,
                                       Phi = theta_baseline, Betas = Beta_hat)
  if (is.null(colnames(theta_baseline))) colnames(theta_baseline) <- colnames(CopyTrack)
  list(R_hat = R_hat, Theta_total = Theta_total, theta_baseline = theta_baseline,
       Beta_hat = Beta_hat, Lambda_hat = Lambda_hat)
}


#' Every metric for one fit
#'
#' RMSE of the signatures, activities, coefficients, intensity and counts, plus
#' either the optimiser's iteration count or the chain's effective sample sizes.
compute_RMSE_paramters <- function(res, data, Lambda_true, Theta_total_true,
                                   R_hat, theta_baseline, Theta_total,
                                   Beta_hat, Lambda_hat, counts_true = NULL) {
  Sens_prec <- Compute_sensitivity_precision(R_hat = R_hat, R_true = data$R)
  MatchedSigs <- match_MutSign(R_true = data$R, R_hat = R_hat)

  rmse_sigs <- sqrt(mean((MatchedSigs$R_hat - MatchedSigs$R_true)^2))
  rmse_theta <- compute_RMSE_Theta(Theta_total_true, Theta_total, MatchedSigs$match)
  rmse_Betas <- compute_RMSE_Betas(data$Betas, Beta_hat, MatchedSigs$match)
  rmse_lambda <- sqrt(mean(rowMeans(Lambda_hat - Lambda_true)^2))
  if (is.null(counts_true)) counts_true <- get_counts_from_data(data)
  rmse_counts <- sqrt(mean(rowMeans(Lambda_hat - counts_true)^2))

  blank <- c(iter = NA_real_, effectiveBetas = NA_real_, effectiveSigs = NA_real_,
             effectiveTheta = NA_real_, effectiveMu = NA_real_,
             effectiveSigma2 = NA_real_, effectiveLogPost = NA_real_,
             effectiveLogLik = NA_real_, effectiveLogPrior = NA_real_)

  if (is.null(res$MCMCchain)) {
    # MAP, or one of the NMF competitors, which report neither.
    sampling_details <- blank
    sampling_details["iter"] <- if (!is.null(res$MAPsolution$it))
      as.numeric(res$MAPsolution$it) else NA_real_
  } else {
    keep_draws <- kept_draw_index(res)
    keep <- colnames(R_hat)
    ess <- function(x) mean(get_PosteriorEffectiveSize(x, keep_draws))
    sampling_details <- c(
      iter = NA_real_,
      effectiveBetas = ess(res$MCMCchain$BETASchain[, , keep, drop = FALSE]),
      effectiveSigs = ess(res$MCMCchain$SIGSchain[, , keep, drop = FALSE]),
      effectiveTheta = ess(res$MCMCchain$THETAchain[, keep, , drop = FALSE]),
      effectiveMu = ess(res$MCMCchain$MUchain[, keep, drop = FALSE]),
      effectiveSigma2 = ess(res$MCMCchain$SIGMA2chain[, keep, drop = FALSE]),
      effectiveLogPost = unname(get_PosteriorEffectiveSize(
        c(res$MCMCchain$logPostchain)[keep_draws], seq_along(keep_draws))),
      effectiveLogLik = unname(get_PosteriorEffectiveSize(
        c(res$MCMCchain$logLikchain)[keep_draws], seq_along(keep_draws))),
      effectiveLogPrior = unname(get_PosteriorEffectiveSize(
        c(res$MCMCchain$logPriorchain)[keep_draws], seq_along(keep_draws))))
  }

  c("Kest" = ncol(R_hat), Sens_prec,
    "rmse_sig" = rmse_sigs, "rmse_theta" = rmse_theta,
    "rmse_Betas" = rmse_Betas, "rmse_lambda" = rmse_lambda,
    "rmse_counts" = rmse_counts,
    "time" = fit_runtime_mins(res),
    sampling_details)
}
