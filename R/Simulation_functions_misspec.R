## Misspecification study: data generation, model fitting and scoring.

################################################################################
# Part 0 - Channel opportunity
################################################################################

if (!exists("get_PosteriorEffectiveSize", mode = "function")) {
  source(file.path(R_DIR, "Utils_functions.R"))
}

compute_mutation_opportunities_hg19 <- function(
    channels = rownames(SignaturePPF::COSMIC_v3.4_SBS96_GRCh37),
    genome = "BSgenome.Hsapiens.UCSC.hg19",
    cache_file = PATH_OPPORTUNITY,
    normalize = TRUE,
    force = FALSE) {
  if (!force && !is.null(cache_file) && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  bs <- BSgenome::getBSgenome(genome)
  # The context of channel "X[Y>Z]W" is the trinucleotide "XYW".
  ctx <- paste0(substr(channels, 1, 1), substr(channels, 3, 3), substr(channels, 7, 7))

  chroms <- paste0("chr", c(1:22, "X", "Y"))
  tri <- NULL
  for (ch in chroms) {
    f <- Biostrings::trinucleotideFrequency(bs[[ch]])
    tri <- if (is.null(tri)) f else tri + f
  }
  # Double-stranded (pyrimidine-centred) count: a context and its reverse
  # complement label the same physical site on opposite strands.
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(names(tri))))
  tri_ds <- tri + tri[match(rc, names(tri))]

  opp <- tri_ds[ctx]
  names(opp) <- channels
  if (normalize) opp <- opp / sum(opp)

  if (!is.null(cache_file)) {
    create_directory(dirname(cache_file))
    saveRDS(opp, cache_file)
  }
  opp
}


#' Draw a channel opportunity multiplier around a given mean
#'
#' Returns a vector summing to `length(opp)`, i.e. mean multiplier one, so the
#' total mutation count stays comparable across scenarios.
sample_opportunity_dirichlet <- function(opp, concentration = 1e5) {
  I <- length(opp)
  stopifnot(concentration > 0, all(opp >= 0), sum(opp) > 0)
  p <- opp / sum(opp)
  draw <- c(LaplacesDemon::rdirichlet(n = 1, alpha = concentration * p))
  opportunity <- I * draw
  names(opportunity) <- names(opp)
  opportunity
}


################################################################################
# Part 1 - The misspecification knobs
################################################################################

#' Configuration of one misspecification scenario
#'
#' Every knob defaults to "no violation", so `misspec_config()` is the correctly
#' specified data-generating process.
misspec_config <- function(epigenome_noise_sd = 0,        # sd of the N(0, sd^2) added to x(t) per patient
                           cn_noise_sd = 0,               # sd of the noise on the OBSERVED copy number
                           cn_min = 0.01,                 # floor for observed / loss copy numbers
                           cn_allow_loss = TRUE,          # deletion states (CN < 2) on by default
                           cn_prob_loss = 0.2,            # probability a segment is a loss
                           cn_loss_states = c(1, 0.5, 0.01),
                           overdispersion_size = Inf,     # NegBin size; Inf => Poisson
                           opportunity = NULL,            # length-96 opportunity vector; NULL => flat
                           opportunity_concentration = 1e5,
                           hotspot_frac = 0,              # fraction of patients carrying hotspots
                           hotspot_n = 5,                 # hotspot tiles per carrier x signature
                           hotspot_width = 0,             # half-width in tiles
                           hotspot_mu = 30,               # expected EXTRA mutations per hotspot tile
                           hotspot_signatures = NULL) {   # NULL => all signatures
  stopifnot(epigenome_noise_sd >= 0, cn_noise_sd >= 0,
            cn_min > 0, cn_prob_loss >= 0, cn_prob_loss <= 1)
  stopifnot(opportunity_concentration > 0)
  stopifnot(is.null(opportunity) || length(opportunity) == 96)
  stopifnot(all(cn_loss_states > 0))
  stopifnot(is.infinite(overdispersion_size) || overdispersion_size > 0)
  stopifnot(hotspot_frac >= 0, hotspot_frac <= 1, hotspot_n >= 0,
            hotspot_width >= 0, hotspot_mu >= 0)
  stopifnot(is.null(hotspot_signatures) || is.character(hotspot_signatures))
  list(epigenome_noise_sd = epigenome_noise_sd,
       cn_noise_sd = cn_noise_sd,
       cn_min = cn_min,
       cn_allow_loss = cn_allow_loss,
       cn_prob_loss = cn_prob_loss,
       cn_loss_states = cn_loss_states,
       overdispersion_size = overdispersion_size,
       opportunity = opportunity,
       opportunity_concentration = opportunity_concentration,
       hotspot_frac = hotspot_frac,
       hotspot_n = hotspot_n,
       hotspot_width = hotspot_width,
       hotspot_mu = hotspot_mu,
       hotspot_signatures = hotspot_signatures)
}


#' Copy number that can also lose material (CN < 2)
#'
#' Returns the TRUE copy number the data is generated from.
generate_CopyTrack_misspec <- function(J = 100,
                                       mu_copy = 1,
                                       size_copy = 10,
                                       length_genome = 2e6,
                                       tilewidth = 100,
                                       mean_segments = 20,
                                       allow_loss = TRUE,
                                       prob_loss = 0.05,
                                       loss_states = c(1, 0.5, 0.01)) {
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
      if (allow_loss && runif(1) < prob_loss) {
        cn <- sample(loss_states, size = 1)
      } else {
        cn <- 2 + rnbinom(1, mu = mu_copy, size = size_copy)
      }
      profile[s:e] <- cn
    }
    cnmat[, j] <- profile
  }
  colnames(cnmat) <- paste0("Patient_", sprintf("%02d", 1:J))
  GenomicRanges::mcols(grX) <- as.data.frame(cnmat)
  grX
}


#' Corrupt a true copy-number matrix, to mimic imperfect CN calling
add_CopyNumber_noise <- function(CopyTrack_true, cn_noise_sd, cn_min = 0.01) {
  if (cn_noise_sd <= 0) return(CopyTrack_true)
  noise <- matrix(rnorm(length(CopyTrack_true), mean = 0, sd = cn_noise_sd),
                  nrow = nrow(CopyTrack_true), ncol = ncol(CopyTrack_true))
  CopyTrack_obs <- CopyTrack_true + noise
  CopyTrack_obs[CopyTrack_obs < cn_min] <- cn_min
  CopyTrack_obs
}


#' Draw the hypermutation hotspots as (patient, signature) records
#'
#' Positions are independent of x(t). `mu` is additive and decoupled from the
#' sparse background, so it sets the height of the spike directly.
draw_hotspot_records <- function(carriers, hot_cols, sig_names, n_tiles,
                                 n_hot = 5, width = 0, mu = 30) {
  empty <- data.frame(patient = integer(0), sample = character(0),
                      signature = character(0), col = integer(0),
                      center = integer(0), start = integer(0), end = integer(0),
                      mu = numeric(0), stringsAsFactors = FALSE)
  if (n_hot <= 0 || mu <= 0 || length(hot_cols) == 0 || !any(carriers)) return(empty)
  recs <- list(); idx <- 1L
  for (j in which(carriers)) {
    for (k in hot_cols) {
      centers <- sample.int(n_tiles, size = min(n_hot, n_tiles))
      recs[[idx]] <- data.frame(
        patient   = j,
        sample    = paste0("Patient_", sprintf("%02d", j)),
        signature = sig_names[k],
        col       = k,
        center    = centers,
        start     = pmax(1L, centers - width),
        end       = pmin(n_tiles, centers + width),
        mu        = mu,
        stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  do.call(rbind, recs)
}


#' True generative attribution probability of every mutation of one patient
#'
#' Proportional to R[i,k] * Theta[k,j] * ExpBetas_j[t,k]. The copy number and
#' tile width cancel across k for the smooth part but are kept anyway, so the
#' ADDITIVE hotspot term is on the same expected-count scale.
compute_true_probs <- function(mut, j, R, Theta, ExpBetas_j, CopyTrack_true,
                               tilewidth, hot_records = NULL) {
  ch_names <- rownames(R)
  tiles <- (mut$pos - 1L) %/% tilewidth + 1L
  chan_idx <- match(mut$channel, ch_names)
  base <- R[chan_idx, , drop = FALSE] * ExpBetas_j[tiles, , drop = FALSE]
  base <- sweep(base, 2, Theta[, j], "*")
  base <- base * (0.5 * tilewidth * CopyTrack_true[tiles, j])
  if (!is.null(hot_records) && nrow(hot_records) > 0) {
    for (r in seq_len(nrow(hot_records))) {
      k0 <- hot_records$col[r]
      inwin <- tiles >= hot_records$start[r] & tiles <= hot_records$end[r]
      if (any(inwin)) base[inwin, k0] <- base[inwin, k0] + hot_records$mu[r] * R[chan_idx[inwin], k0]
    }
  }
  rs <- rowSums(base); rs[!is.finite(rs) | rs <= 0] <- 1
  base / rs
}


#' All mutations of one patient, as a superposition of per-signature processes
sample_mutations_patient_misspec <- function(j, R, Theta, ExpBetas_j, CopyTrack_true,
                                             tilewidth = 100,
                                             opportunity = 1,
                                             nb_size = Inf,
                                             hot_records = NULL) {
  CumSums <- apply(ExpBetas_j, 2, function(x) cumsum(0.5 * tilewidth * CopyTrack_true[, j] * x))
  Tot <- CumSums[nrow(CumSums), ]
  I <- nrow(R)
  K <- ncol(R)
  ch_names <- rownames(R)
  sig_names <- colnames(R)
  pieces <- list()
  idx <- 1L
  for (i in 1:I) {
    opp_i <- if (length(opportunity) == 1) opportunity else opportunity[i]
    for (k in 1:K) {
      lam_tot <- opp_i * R[i, k] * Theta[k, j] * Tot[k]
      if (!is.finite(lam_tot) || lam_tot <= 0) next
      if (is.null(nb_size) || is.infinite(nb_size)) {
        N <- rpois(1, lam_tot)
      } else {
        N <- rnbinom(1, mu = lam_tot, size = nb_size)
      }
      if (N > 0) {
        u <- sort(runif(N))
        ts <- findInterval(u, CumSums[, k] / Tot[k]) * tilewidth +
          sample.int(tilewidth - 1, size = N, replace = TRUE)
        pieces[[idx]] <- data.frame(pos = ts,
                                    channel = ch_names[i],
                                    signature = sig_names[k],
                                    stringsAsFactors = FALSE)
        idx <- idx + 1L
      }
    }
  }
  if (!is.null(hot_records) && nrow(hot_records) > 0) {
    for (r in seq_len(nrow(hot_records))) {
      k <- hot_records$col[r]
      Nx <- rpois(1, hot_records$mu[r])
      if (Nx <= 0) next
      chan <- sample.int(I, size = Nx, replace = TRUE, prob = R[, k])
      tiles <- if (hot_records$start[r] == hot_records$end[r]) rep(hot_records$start[r], Nx)
               else sample(hot_records$start[r]:hot_records$end[r], size = Nx, replace = TRUE)
      pos <- (tiles - 1L) * tilewidth + sample.int(tilewidth - 1, size = Nx, replace = TRUE)
      pieces[[idx]] <- data.frame(pos = pos,
                                  channel = ch_names[chan],
                                  signature = sig_names[k],
                                  stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
  sig_cols <- paste0("ptrue_", sig_names)
  if (length(pieces) == 0) {
    empty <- data.frame(pos = integer(0), channel = character(0),
                        signature = character(0), stringsAsFactors = FALSE)
    empty[sig_cols] <- rep(list(numeric(0)), length(sig_cols))
    return(empty)
  }
  mut <- do.call(rbind, pieces)
  Pt <- compute_true_probs(mut, j, R, Theta, ExpBetas_j, CopyTrack_true, tilewidth, hot_records)
  colnames(Pt) <- sig_cols
  cbind(mut, as.data.frame(Pt))
}


#' Sample every patient
sample_dataset_misspec <- function(R, Theta, Xcovs, Betas, ExpBetas_shared,
                                   CopyTrack_true, tilewidth,
                                   opportunity, nb_size,
                                   epigenome_noise_sd,
                                   hotspot_frac = 0, hotspot_n = 5,
                                   hotspot_width = 0, hotspot_mu = 30,
                                   hotspot_signatures = NULL,
                                   ncores = 1,
                                   verbose = FALSE) {
  J <- ncol(Theta)
  n_tiles <- nrow(CopyTrack_true)

  hot_cols <- if (is.null(hotspot_signatures)) seq_len(ncol(R)) else which(colnames(R) %in% hotspot_signatures)

  carriers <- logical(J)
  n_carriers <- round(hotspot_frac * J)
  if (n_carriers > 0) carriers[sample.int(J, size = n_carriers)] <- TRUE
  hotspot_records <- draw_hotspot_records(carriers = carriers, hot_cols = hot_cols,
                                          sig_names = colnames(R), n_tiles = n_tiles,
                                          n_hot = hotspot_n, width = hotspot_width,
                                          mu = hotspot_mu)

  if (ncores > 1) doParallel::registerDoParallel(ncores) else foreach::registerDoSEQ()
  j <- NULL   # keeps R CMD check and linters quiet about the foreach variable
  Mutations <- foreach::foreach(j = 1:J, .combine = "rbind") %dopar% {
    if (verbose) message("Simulating Patient_", sprintf("%02d", j))
    if (epigenome_noise_sd > 0) {
      Xj <- Xcovs + matrix(rnorm(length(Xcovs), mean = 0, sd = epigenome_noise_sd),
                           nrow = nrow(Xcovs), ncol = ncol(Xcovs))
      Xj <- scale(Xj)
      ExpBetas_j <- exp(Xj %*% Betas)
    } else {
      ExpBetas_j <- ExpBetas_shared
    }
    df <- sample_mutations_patient_misspec(
      j = j, R = R, Theta = Theta,
      ExpBetas_j = ExpBetas_j,
      CopyTrack_true = CopyTrack_true,
      tilewidth = tilewidth,
      opportunity = opportunity,
      nb_size = nb_size,
      hot_records = hotspot_records[hotspot_records$patient == j, , drop = FALSE])
    df$sample <- paste0("Patient_", sprintf("%02d", j))
    df
  }
  Mutations$sample <- factor(Mutations$sample, levels = unique(Mutations$sample))
  Mutations$channel <- as.factor(Mutations$channel)
  attr(Mutations, "hotspot_records") <- hotspot_records
  attr(Mutations, "hotspot_carriers") <- which(carriers)
  Mutations
}


#' One misspecification dataset
generate_MutationData_misspec <- function(J = 100,
                                          cosmic_sigs = c("SBS1", "SBS2", "SBS3",
                                                          "SBS5", "SBS8", "SBS13",
                                                          "SBS18", "SBS44"),
                                          K_new = 0,
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
                                          misspec = misspec_config(),
                                          ncores = 1) {

  # Step 1 - covariance of the covariate track
  if (corr == "indep") {
    Sigma <- diag(p_all)
  } else if (corr == "onion") {
    corrmat <- clusterGeneration::genPositiveDefMat(p_all, covMethod = "onion")$Sigma
    Sigma <- cov2cor(corrmat)
  }

  # Step 2 - shared reference epigenome x(t)
  grX <- generate_SignalTrack(length_genome = length_genome,
                              tilewidth = tilewidth, rho = rho, p = p_all,
                              Sigma = Sigma)

  # Step 3 - true per-patient copy number (possibly with deletions)
  gr_CopyTrack_true <- generate_CopyTrack_misspec(J = J, mu_copy = mu_copy,
                                                  size_copy = size_copy,
                                                  length_genome = length_genome,
                                                  tilewidth = tilewidth,
                                                  allow_loss = misspec$cn_allow_loss,
                                                  prob_loss = misspec$cn_prob_loss,
                                                  loss_states = misspec$cn_loss_states)

  # Step 4 - model parameters
  pars_true <- generate_Parameters(cosmic_sigs = cosmic_sigs,
                                   K_new = K_new, J = J,
                                   a = a,
                                   p_all = p_all,
                                   p_to_use = p_to_use,
                                   theta = theta,
                                   sd_beta = sd_beta)

  # Step 5 - channel opportunity (mean 1, so counts stay comparable)
  I <- nrow(pars_true$R)
  if (is.null(misspec$opportunity)) {
    opportunity <- rep(1, I)
    names(opportunity) <- rownames(pars_true$R)
  } else {
    stopifnot(length(misspec$opportunity) == I)
    opp <- misspec$opportunity
    names(opp) <- rownames(pars_true$R)
    opportunity <- sample_opportunity_dirichlet(
      opp = opp, concentration = misspec$opportunity_concentration)
  }

  # Step 6 - rescale tracks and hyperparameters
  Xcovs <- as.matrix(GenomicRanges::mcols(grX)[, -1])   # drop bin_weight
  bin_weight <- grX$bin_weight
  ExpBetas <- exp(Xcovs %*% pars_true$Betas)
  CopyTrack_true <- as.matrix(GenomicRanges::mcols(gr_CopyTrack_true))
  colnames(pars_true$Theta) <- colnames(CopyTrack_true)
  Theta_scaled <- pars_true$Theta / sum(bin_weight)

  # Step 7 - simulate from the TRUE copy number and the (perturbed) epigenome
  nb_size <- misspec$overdispersion_size
  Mutations <- sample_dataset_misspec(R = pars_true$R,
                                      Theta = Theta_scaled,
                                      Xcovs = Xcovs,
                                      Betas = pars_true$Betas,
                                      ExpBetas_shared = ExpBetas,
                                      CopyTrack_true = CopyTrack_true,
                                      tilewidth = tilewidth,
                                      opportunity = opportunity,
                                      nb_size = nb_size,
                                      epigenome_noise_sd = misspec$epigenome_noise_sd,
                                      hotspot_frac = misspec$hotspot_frac,
                                      hotspot_n = misspec$hotspot_n,
                                      hotspot_width = misspec$hotspot_width,
                                      hotspot_mu = misspec$hotspot_mu,
                                      hotspot_signatures = misspec$hotspot_signatures,
                                      ncores = ncores)
  hotspot_records <- attr(Mutations, "hotspot_records")
  hotspot_carriers <- attr(Mutations, "hotspot_carriers")
  if (nrow(hotspot_records) > 0) {
    hotspot_records$start_bp <- (hotspot_records$start - 1L) * tilewidth + 1L
    hotspot_records$end_bp <- hotspot_records$end * tilewidth
  }

  ptrue_cols <- grep("^ptrue_", colnames(Mutations), value = TRUE)
  p_true <- as.matrix(Mutations[, ptrue_cols, drop = FALSE])
  colnames(p_true) <- sub("^ptrue_", "", ptrue_cols)
  p_true <- p_true[, colnames(pars_true$R), drop = FALSE]

  # Step 8 - the mutation GRanges, carrying the SHARED covariate values
  gr_Mutations <- GenomicRanges::GRanges(
    seqnames = "chrsim",
    IRanges::IRanges(start = Mutations$pos, end = Mutations$pos),
    sample = Mutations$sample,
    channel = Mutations$channel)
  overlaps <- GenomicRanges::findOverlaps(gr_Mutations, grX)
  GenomicRanges::mcols(gr_Mutations) <- cbind(
    GenomicRanges::mcols(gr_Mutations),
    GenomicRanges::mcols(grX[S4Vectors::subjectHits(overlaps)]))
  gr_Mutations$bin_weight <- NULL
  # Kept OUT of mcols, so the fitting call sees exactly sample, channel and the
  # covariates - which is also what SignaturePPF_validate() requires.
  signature_true <- as.character(Mutations$signature)

  # Step 9 - the observed (noisy) copy number the model will actually use
  CopyTrack_obs <- add_CopyNumber_noise(CopyTrack_true,
                                        cn_noise_sd = misspec$cn_noise_sd,
                                        cn_min = misspec$cn_min)
  gr_CopyTrack_obs <- gr_CopyTrack_true
  GenomicRanges::mcols(gr_CopyTrack_obs) <- as.data.frame(CopyTrack_obs)
  Betas <- pars_true$Betas
  rownames(Betas) <- colnames(Xcovs)
  colnames(Betas) <- colnames(pars_true$R)

  list("gr_Mutations" = gr_Mutations,
       "signature_true" = signature_true,
       "p_true" = p_true,
       "covariate_used" = colnames(GenomicRanges::mcols(grX)[, -1])[1:p_to_use],
       "R" = pars_true$R,
       "Betas" = Betas,
       "Theta" = Theta_scaled,
       "opportunity" = opportunity,
       "gr_SignalTrack" = grX,
       "gr_CopyTrack" = gr_CopyTrack_obs,        # observed (model input)
       "gr_CopyTrack_true" = gr_CopyTrack_true,  # truth
       "CopyTrack_true" = CopyTrack_true,
       "hotspot_records" = hotspot_records,
       "hotspot_carriers" = hotspot_carriers,
       "misspec" = misspec,
       "simulation_parameters" = list(J = J,
                                      cosmic_sigs = cosmic_sigs,
                                      K_new = K_new,
                                      theta = theta,
                                      a = a,
                                      mu_copy = mu_copy,
                                      size_copy = size_copy,
                                      length_genome = length_genome,
                                      tilewidth = tilewidth,
                                      rho = rho,
                                      p_all = p_all,
                                      p_to_use = p_to_use,
                                      corr = corr))
}


################################################################################
# Part 2 - Fitting
################################################################################

#' Train / test split of the genome
split_bins_misspec <- function(data, n_test_bins = 4000) {
  n_bins <- length(data$gr_SignalTrack)
  n_train <- n_bins - n_test_bins
  if (n_train <= 0) stop("n_test_bins must be smaller than the number of tiles")
  bin_mut <- S4Vectors::subjectHits(
    GenomicRanges::findOverlaps(data$gr_Mutations, data$gr_SignalTrack))
  list(n_bins = n_bins,
       n_train = n_train,
       train_bins = seq_len(n_train),
       test_bins = seq.int(n_train + 1L, n_bins),
       bin_mut = bin_mut,
       train_mut = which(bin_mut <= n_train),
       test_mut = which(bin_mut > n_train))
}


#' The training tracks and mutations, in SignaturePPF form
#'
#' The observed (noisy) copy number is what the models get - exactly what an
#' analyst would have. Copy number carries the tile width and the factor of one
#' half, so the intensity is per tile rather than per base.
ppf_training_data <- function(data, n_test_bins = 4000) {
  tilewidth <- data$simulation_parameters$tilewidth
  sp <- split_bins_misspec(data, n_test_bins = n_test_bins)
  list(
    gr_Mutations = data$gr_Mutations[sp$train_mut],
    SignalTrack = as.matrix(GenomicRanges::mcols(data$gr_SignalTrack))[sp$train_bins, -1],
    CopyTrack = tilewidth *
      as.matrix(GenomicRanges::mcols(data$gr_CopyTrack))[sp$train_bins, ] / 2)
}


#' DE NOVO fits: the two competitors plus SignaturePPF, signatures ESTIMATED
run_models_misspec <- function(out_dir,
                               K = 15,
                               n_test_bins = 4000,
                               run_CompNMF = TRUE,
                               run_SignatureAnalyzer = TRUE,
                               run_MAP = TRUE,
                               run_MCMC = FALSE,
                               overwrite = FALSE,
                               seed = 1L,
                               a = 1.01, alpha = 1.01, epsilon = 0.001,
                               c0 = 100, d0 = 1,
                               tol = 1e-6,
                               maxiter = 4000,
                               nsamples = 3000,
                               burnin = 1500,
                               thin = 1,
                               sa_K0 = 15,
                               sa_nrun = 5,
                               sa_niter = 1e5) {
  data <- readRDS(file.path(out_dir, "data.rds.gzip"))
  train <- ppf_training_data(data, n_test_bins = n_test_bins)
  MutMatrix <- SignaturePPF::getTotalMutations(train$gr_Mutations)

  controls <- SignaturePPF::SignaturePPF_control(
    maxiter = maxiter, tol = tol, nsamples = nsamples, burnin = burnin,
    thin = thin, print_every = 0L)
  prior <- SignaturePPF::SignaturePPF_prior(a = a, alpha = alpha,
                                            epsilon = epsilon, c0 = c0, d0 = d0)

  todo <- function(file) overwrite || !file.exists(file.path(out_dir, file))

  # ---- Competitor 1: compressive NMF (no covariates, no copy number) ----
  if (run_CompNMF && todo("output_CompNMF.rds.gzip")) {
    t0 <- Sys.time()
    outCompNMF <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = K, a = a,
                                                     alpha = alpha,
                                                     epsilon = epsilon, tol = 1e-7)
    outCompNMF$time <- Sys.time() - t0
    saveRDS(outCompNMF, file.path(out_dir, "output_CompNMF.rds.gzip"), compress = "gzip")
  }

  # ---- Competitor 2: SignatureAnalyzer / BayesNMF ----
  if (run_SignatureAnalyzer && todo("output_SignatureAnalyzer.rds.gzip")) {
    t0 <- Sys.time()
    outSA <- sigminer::sig_auto_extract(nmf_matrix = t(MutMatrix), K0 = sa_K0,
                                        nrun = sa_nrun, niter = sa_niter,
                                        cores = 1, destdir = tempfile("SA_"))
    outSA$time <- Sys.time() - t0
    saveRDS(outSA, file.path(out_dir, "output_SignatureAnalyzer.rds.gzip"), compress = "gzip")
  }

  # ---- Proposed: SignaturePPF, MAP ----
  if (run_MAP && todo("output_map_FullModel.rds.gzip")) {
    out_map <- SignaturePPF::SignaturePPF(
      train, sigs = NULL, K = K, method = "map",
      prior = prior, controls = controls, seed = seed, verbose = FALSE,
      prune_solution = FALSE)   # scored by this file's own rule; see Simulation_main.R
    saveRDS(out_map, file.path(out_dir, "output_map_FullModel.rds.gzip"), compress = "gzip")
  }

  # ---- Proposed: SignaturePPF, MCMC started at the MAP ----
  if (run_MCMC && todo("output_mcmc_FullModel.rds.gzip")) {
    out_map <- readRDS(file.path(out_dir, "output_map_FullModel.rds.gzip"))
    init <- SignaturePPF::SignaturePPF_init(
      R_start = out_map$Signatures,
      Theta_start = out_map$Thetas,        # activity scale, as the package expects
      Betas_start = out_map$Betas,
      Mu_start = out_map$Mu,
      Sigma2_start = out_map$Sigma2)
    out_mcmc <- SignaturePPF::SignaturePPF(
      train, sigs = NULL, K = K, method = "mcmc",
      prior = prior, controls = controls, init = init,
      init_mcmc_from_map = FALSE, seed = seed, verbose = FALSE,
      prune_solution = FALSE)
    saveRDS(out_mcmc, file.path(out_dir, "output_mcmc_FullModel.rds.gzip"), compress = "gzip")
  }
  invisible(NULL)
}


#' FIXED-signature fits, for the attribution / calibration use case
run_models_misspec_fixed <- function(out_dir,
                                     distractors = c("SBS6", "SBS20", "SBS26",
                                                     "SBS40a", "SBS30"),
                                     n_test_bins = 4000,
                                     run_CompNMF = TRUE,
                                     run_MAP = TRUE,
                                     run_MCMC = FALSE,
                                     overwrite = FALSE,
                                     seed = 1L,
                                     a = 1.01, alpha = 1.01, epsilon = 0.001,
                                     c0 = 100, d0 = 1,
                                     tol = 1e-6,
                                     maxiter = 4000,
                                     nsamples = 3000,
                                     burnin = 1500,
                                     thin = 1) {
  data <- readRDS(file.path(out_dir, "data.rds.gzip"))
  train <- ppf_training_data(data, n_test_bins = n_test_bins)

  controls <- SignaturePPF::SignaturePPF_control(
    maxiter = maxiter, tol = tol, nsamples = nsamples, burnin = burnin,
    thin = thin, print_every = 0L)
  prior <- SignaturePPF::SignaturePPF_prior(a = a, alpha = alpha,
                                            epsilon = epsilon, c0 = c0, d0 = d0)

  catalog <- unique(c(colnames(data$R), distractors))
  R_fixed <- as.matrix(SignaturePPF::COSMIC_v3.4_SBS96_GRCh37[, catalog])

  todo <- function(file) overwrite || !file.exists(file.path(out_dir, file))

  # ---- Fixed-signature CompNMF: signatures pinned through a concentrated prior.
  #      The position-independent benchmark for both counts and calibration.
  if (run_CompNMF && todo("output_CompNMF_Fixed.rds.gzip")) {
    MutMatrix <- SignaturePPF::getTotalMutations(train$gr_Mutations)
    t0 <- Sys.time()
    outCompNMF <- CompressiveNMF::CompressiveNMF_map(MutMatrix, K = 0,
                                                     S = 1e7 * R_fixed + 1)
    outCompNMF$time <- Sys.time() - t0
    saveRDS(outCompNMF, file.path(out_dir, "output_CompNMF_Fixed.rds.gzip"), compress = "gzip")
  }

  if (run_MAP && todo("output_map_Fixed.rds.gzip")) {
    out_map <- SignaturePPF::SignaturePPF(
      train, sigs = R_fixed, sigs_fixed = TRUE, method = "map",
      prior = prior, controls = controls, seed = seed, verbose = FALSE,
      prune_solution = FALSE)   # a refit reports every reference, parked or not
    saveRDS(out_map, file.path(out_dir, "output_map_Fixed.rds.gzip"), compress = "gzip")
  }

  if (run_MCMC && todo("output_mcmc_Fixed.rds.gzip")) {
    out_map <- readRDS(file.path(out_dir, "output_map_Fixed.rds.gzip"))
    # No R_start here: in a refit the signatures come from `sigs`, and passing
    # both is rejected by the package rather than silently ignored.
    init <- SignaturePPF::SignaturePPF_init(
      Theta_start = out_map$Thetas,
      Betas_start = out_map$Betas,
      Mu_start = out_map$Mu,
      Sigma2_start = out_map$Sigma2)
    out_mcmc <- SignaturePPF::SignaturePPF(
      train, sigs = R_fixed, sigs_fixed = TRUE, method = "mcmc",
      prior = prior, controls = controls, init = init,
      init_mcmc_from_map = FALSE, seed = seed, verbose = FALSE,
      prune_solution = FALSE)
    saveRDS(out_mcmc, file.path(out_dir, "output_mcmc_Fixed.rds.gzip"), compress = "gzip")
  }
  invisible(NULL)
}


################################################################################
# Part 3 - Scoring
################################################################################

#' Runtime of a fit, in minutes
#'
#' SignaturePPF records `$runtime` itself; the competitors get a `$time` set by
#' the fitting functions above.
fit_runtime_mins <- function(res) {
  dt <- if (!is.null(res$runtime)) res$runtime else res$time
  if (is.null(dt)) NA_real_ else as.numeric(dt, units = "mins")
}


#' The optimiser's iteration count, or the chain's effective sample sizes
sampling_details_of <- function(res, keep = NULL) {
  blank <- c(iter = NA_real_, effectiveBetas = NA_real_, effectiveSigs = NA_real_,
             effectiveTheta = NA_real_, effectiveMu = NA_real_,
             effectiveSigma2 = NA_real_, effectiveLogPost = NA_real_,
             effectiveLogLik = NA_real_, effectiveLogPrior = NA_real_)

  if (is.null(res$MCMCchain)) {
    it <- if (!is.null(res$MAPsolution$iter)) res$MAPsolution$iter
          else if (!is.null(res$mapOutput$iter)) res$mapOutput$iter
          else NA_real_
    blank["iter"] <- as.numeric(it)
    return(blank)
  }

  keep_draws <- kept_draw_index(res)
  chain_sigs <- dimnames(res$MCMCchain$MUchain)[[2]]
  if (is.null(keep)) keep <- chain_sigs else keep <- intersect(keep, chain_sigs)
  if (!length(keep)) keep <- chain_sigs
  ess <- function(x) mean(get_PosteriorEffectiveSize(x, keep_draws), na.rm = TRUE)
  scalar_ess <- function(x) unname(get_PosteriorEffectiveSize(
    c(x)[keep_draws], seq_along(keep_draws)))

  c(iter = NA_real_,
    effectiveBetas = ess(res$MCMCchain$BETASchain[, , keep, drop = FALSE]),
    effectiveSigs = ess(res$MCMCchain$SIGSchain[, , keep, drop = FALSE]),
    effectiveTheta = ess(res$MCMCchain$THETAchain[, keep, , drop = FALSE]),
    effectiveMu = ess(res$MCMCchain$MUchain[, keep, drop = FALSE]),
    effectiveSigma2 = ess(res$MCMCchain$SIGMA2chain[, keep, drop = FALSE]),
    effectiveLogPost = scalar_ess(res$MCMCchain$logPostchain),
    effectiveLogLik = scalar_ess(res$MCMCchain$logLikchain),
    effectiveLogPrior = scalar_ess(res$MCMCchain$logPriorchain))
}


#' A common estimate structure for any of the three model types
extract_estimates_misspec <- function(res, model_type, SignalTrack, CopyTrack) {
  flat <- matrix(rep(1 / 96, 96))
  if (model_type == "PPF") {
    cutoff <- 5 * res$prior$a * res$prior$epsilon
    keep <- (res$Mu > 5 * cutoff) &
      (c(sigminer::cosine(res$Signatures, flat)) < 0.975)
    if (!any(keep)) keep <- rep(TRUE, length(keep))
    R_hat <- res$Signatures[, keep, drop = FALSE]
    Beta_hat <- res$Betas[, keep, drop = FALSE]
    # $Baseline, NOT $Thetas: under the activity prior the latter is the total.
    theta_baseline <- res$Baseline[keep, , drop = FALSE]
    ExpXB <- exp(SignalTrack[, rownames(Beta_hat), drop = FALSE] %*% Beta_hat)
    Theta_total <- theta_baseline * crossprod(ExpXB, CopyTrack)
    Lambda_hat <- reconstruct_lambda_raw(
      SignalTrack = SignalTrack, CopyTrack = CopyTrack,
      Phi = theta_baseline, Betas = Beta_hat)
  } else {
    if (model_type == "CompNMF") {
      W <- res$Signatures; H <- res$Theta
      keep <- (res$Mu > 0) & (c(sigminer::cosine(W, flat)) < 0.975)
    } else if (model_type == "SignatureAnalyzer") {
      W <- res$Signature.norm; H <- res$Exposure
      keep <- c(sigminer::cosine(W, flat)) < 0.975
    } else {
      stop("Unknown model_type: ", model_type)
    }
    if (!any(keep)) keep <- rep(TRUE, length(keep))
    R_hat <- W[, keep, drop = FALSE]
    Theta_total <- H[keep, , drop = FALSE]
    # Per-copy baseline, so the intensity re-applies the copy number
    denom <- t(colSums(CopyTrack))[rep(1, nrow(Theta_total)), , drop = FALSE]
    theta_baseline <- Theta_total / denom
    Beta_hat <- matrix(0, nrow = ncol(SignalTrack), ncol = nrow(theta_baseline),
                       dimnames = list(colnames(SignalTrack), rownames(theta_baseline)))
    Lambda_hat <- reconstruct_lambda_raw(
      SignalTrack = SignalTrack, CopyTrack = CopyTrack,
      Phi = theta_baseline, Betas = Beta_hat)
  }
  if (is.null(colnames(theta_baseline))) colnames(theta_baseline) <- colnames(CopyTrack)
  list(R_hat = R_hat, theta_baseline = theta_baseline, Beta_hat = Beta_hat,
       Theta_total = Theta_total, Lambda_hat = Lambda_hat)
}


#' Posterior probability that each mutation came from each estimated signature
compute_assignment_probs <- function(R_hat, theta_baseline, Beta_hat, X_mut,
                                     channel_id, sample_id) {
  n <- length(channel_id); K <- ncol(R_hat)
  num <- R_hat[channel_id, , drop = FALSE] *
    t(theta_baseline[, sample_id, drop = FALSE])
  if (!(is.null(Beta_hat) || is.null(rownames(Beta_hat)) || all(Beta_hat == 0))) {
    num <- num * exp(as.matrix(X_mut[, rownames(Beta_hat), drop = FALSE]) %*% Beta_hat)
  }
  rs <- rowSums(num)
  bad <- !is.finite(rs) | rs <= 0
  rs[bad] <- 1
  P <- num / rs
  if (any(bad)) P[bad, ] <- 1 / K
  colnames(P) <- colnames(R_hat)
  P
}


#' Label each estimated signature with the true one it best matches
match_hat_to_true_cols <- function(R_true, R_hat, cos_cutoff = 0.8) {
  k_true <- ncol(R_true); k_hat <- ncol(R_hat)
  ms <- match_MutSign(R_true, R_hat)
  lab <- rep(NA_character_, k_hat)
  for (i in seq_len(k_true)) {
    h <- ms$match[i]
    if (h <= k_hat && sigminer::cosine(R_true[, i], R_hat[, h]) >= cos_cutoff) {
      lab[h] <- colnames(R_true)[i]
    }
  }
  lab
}


#' Expected / maximum calibration error and the reliability curve
compute_ECE <- function(prob, correct, nbins = 10) {
  edges <- seq(0, 1, length.out = nbins + 1)
  bin <- findInterval(prob, edges, rightmost.closed = TRUE, all.inside = TRUE)
  n <- length(prob); ece <- 0; mce <- 0
  rel <- data.frame()
  for (b in seq_len(nbins)) {
    idx <- which(bin == b)
    if (length(idx) == 0) next
    conf_b <- mean(prob[idx]); acc_b <- mean(correct[idx]); w <- length(idx) / n
    gap <- abs(acc_b - conf_b)
    ece <- ece + w * gap; mce <- max(mce, gap)
    rel <- rbind(rel, data.frame(bin = b, mid = (edges[b] + edges[b + 1]) / 2,
                                 mean_conf = conf_b, accuracy = acc_b, n = length(idx)))
  }
  list(ece = ece, mce = mce, reliability = rel)
}


#' Mutation-level attribution accuracy and probability calibration
evaluate_attribution_calibration <- function(P, lab, signature_true, true_sigs,
                                             nbins = 10, p_oracle = NULL) {
  n <- nrow(P)
  amax <- max.col(P, ties.method = "first")
  pred_label <- lab[amax]
  correct <- !is.na(pred_label) & (pred_label == signature_true)
  accuracy <- mean(correct)
  conf <- P[cbind(seq_len(n), amax)]
  cal <- compute_ECE(conf, correct, nbins = nbins)

  # Macro-averaged over the true signatures, so a rare signature counts as much
  # as an abundant one.
  prec_s <- rec_s <- f1_s <- rep(NA_real_, length(true_sigs))
  for (si in seq_along(true_sigs)) {
    s <- true_sigs[si]
    tp <- sum(!is.na(pred_label) & pred_label == s & signature_true == s)
    pp <- sum(!is.na(pred_label) & pred_label == s)
    ap <- sum(signature_true == s)
    prec_s[si] <- if (pp > 0) tp / pp else 0
    rec_s[si] <- if (ap > 0) tp / ap else NA_real_
    denom <- prec_s[si] + rec_s[si]
    f1_s[si] <- if (!is.na(denom) && denom > 0) 2 * prec_s[si] * rec_s[si] / denom else 0
  }
  precision <- mean(prec_s, na.rm = TRUE)
  sensitivity <- mean(rec_s, na.rm = TRUE)
  f1 <- mean(f1_s, na.rm = TRUE)

  P_model <- matrix(0, n, length(true_sigs), dimnames = list(NULL, true_sigs))
  for (s in true_sigs) {
    cols <- which(lab == s)
    if (length(cols) > 0) P_model[, s] <- rowSums(P[, cols, drop = FALSE])
  }
  Y <- outer(signature_true, true_sigs, FUN = "==") * 1
  brier <- mean(rowSums((P_model - Y)^2))

  # Oracle calibration against the TRUE generative probabilities. This is the
  # cleanest target under misspecification: it does not depend on the single
  # noisy label that was actually drawn.
  oracle_l2 <- NA_real_; oracle_kl <- NA_real_; estim_error <- NA_real_
  if (!is.null(p_oracle)) {
    p_oracle <- p_oracle[, true_sigs, drop = FALSE]
    oracle_l2 <- mean(rowSums((P_model - p_oracle)^2))
    eps <- 1e-12
    oracle_kl <- mean(rowSums(p_oracle * (log(p_oracle + eps) - log(P_model + eps))))

    # Fraction of mutations attributed differently from the attribution the TRUE
    # parameters would make. Unlike accuracy against the realised label this is
    # not capped by the Bayes error of overlapping signatures: a perfectly
    # estimated model scores 0 whatever the irreducible uncertainty.
    oracle_label <- true_sigs[max.col(p_oracle, ties.method = "first")]
    estim_error <- mean(is.na(pred_label) | pred_label != oracle_label)
  }
  list(accuracy = accuracy, precision = precision, sensitivity = sensitivity, f1 = f1,
       brier = brier, ece = cal$ece, mce = cal$mce,
       oracle_l2 = oracle_l2, oracle_kl = oracle_kl, estim_error = estim_error,
       reliability = cal$reliability)
}


#' Which covariate effects does the model declare non-zero?
beta_CI_excludes_zero <- function(fit, level = 0.95, dimnames_ref = NULL) {
  if (is.null(fit$MCMCchain) || is.null(fit$MCMCchain$BETASchain)) return(NULL)
  ci <- SignaturePPF::posterior_CI(fit, what = "Betas", level = level)
  out <- (ci$lowCI > 0) | (ci$highCI < 0)
  if (!is.null(dimnames_ref)) dimnames(out) <- dimnames_ref
  out
}


#' False-positive rate over the truly-zero covariate effects
compute_FP_Betas <- function(Beta_true, Beta_hat, tol = 0.05, selected = NULL) {
  if (is.null(Beta_hat) || is.null(dimnames(Beta_hat))) return(NA_real_)
  sigs <- intersect(colnames(Beta_true), colnames(Beta_hat))
  covs <- intersect(rownames(Beta_true), rownames(Beta_hat))
  if (!length(sigs) || !length(covs)) return(NA_real_)
  bt <- Beta_true[covs, sigs, drop = FALSE]
  bh <- Beta_hat[covs, sigs, drop = FALSE]
  zero <- (bt == 0)
  if (!any(zero)) return(NA_real_)
  decl <- if (is.null(selected)) (abs(bh) > tol) else selected[covs, sigs, drop = FALSE]
  mean(decl[zero])
}


#' Sign concordance over the truly non-zero covariate effects
compute_sign_Betas <- function(Beta_true, Beta_hat, selected = NULL) {
  if (is.null(Beta_hat) || is.null(dimnames(Beta_hat))) return(NA_real_)
  sigs <- intersect(colnames(Beta_true), colnames(Beta_hat))
  covs <- intersect(rownames(Beta_true), rownames(Beta_hat))
  if (!length(sigs) || !length(covs)) return(NA_real_)
  bt <- Beta_true[covs, sigs, drop = FALSE]
  bh <- Beta_hat[covs, sigs, drop = FALSE]
  keep <- (bt != 0)
  if (!is.null(selected)) keep <- keep & selected[covs, sigs, drop = FALSE]
  if (!any(keep)) return(NA_real_)
  mean(sign(bh[keep]) == sign(bt[keep]))
}


#' DE NOVO evaluation: reconstruction of signatures, activities, effects, counts
postProcessOutput_misspec <- function(out_dir, beta_tol = 0.05, n_test_bins = 4000) {
  data <- readRDS(file.path(out_dir, "data.rds.gzip"))
  tilewidth <- data$simulation_parameters$tilewidth
  SignalTrackSim <- as.matrix(GenomicRanges::mcols(data$gr_SignalTrack))[, -1]
  CopyTrack_obs <- tilewidth * as.matrix(GenomicRanges::mcols(data$gr_CopyTrack)) / 2
  CopyTrack_true <- tilewidth * data$CopyTrack_true / 2

  bsp <- split_bins_misspec(data, n_test_bins = n_test_bins)

  Lambda_true <- reconstruct_lambda_raw(
    SignalTrack = SignalTrackSim, CopyTrack = CopyTrack_true,
    Phi = data$Theta, Betas = data$Betas)
  Theta_total_true <- data$Theta * crossprod(
    exp(SignalTrackSim[bsp$train_bins, , drop = FALSE] %*% data$Betas),
    CopyTrack_true[bsp$train_bins, , drop = FALSE])

  sm_mut <- as.character(data$gr_Mutations$sample)
  n_bins <- nrow(SignalTrackSim)
  mat_counts <- unclass(table(factor(bsp$bin_mut, levels = seq_len(n_bins)),
                              factor(sm_mut, levels = colnames(CopyTrack_obs))))
  storage.mode(mat_counts) <- "double"

  lambda_on_bins <- function(est, bins) reconstruct_lambda_raw(
    SignalTrack = SignalTrackSim[bins, , drop = FALSE],
    CopyTrack = CopyTrack_obs[bins, , drop = FALSE],
    Phi = est$theta_baseline, Betas = est$Beta_hat)

  eval_one <- function(file, model_type, model_name) {
    res <- open_rds_file(file.path(out_dir, file))
    if (is.null(res)) return(NULL)
    # Estimates are extracted on the TRAINING tracks, so the position-independent
    # competitors spread their exposures over the fitted region only.
    est <- extract_estimates_misspec(res, model_type,
                                     SignalTrackSim[bsp$train_bins, , drop = FALSE],
                                     CopyTrack_obs[bsp$train_bins, , drop = FALSE])

    sp <- Compute_sensitivity_precision(est$R_hat, data$R)
    ms <- match_MutSign(data$R, est$R_hat)
    rmse_sig <- sqrt(mean((ms$R_hat - ms$R_true)^2))
    cosine_R <- mean(sapply(seq_len(ncol(data$R)), function(s)
      max(apply(est$R_hat, 2, function(h) sigminer::cosine(data$R[, s], h)))))
    rmse_theta <- compute_RMSE_Theta(Theta_total_true, est$Theta_total, ms$match)
    rmse_Betas <- compute_RMSE_Betas(data$Betas, est$Beta_hat, ms$match)

    # Relabel the estimated Beta columns with the true signature they matched,
    # so the truly-zero entries can be scored by name.
    Bhat <- est$Beta_hat
    lab_cols <- rep(NA_character_, ncol(Bhat))
    for (s in seq_len(ncol(data$R))) {
      h <- ms$match[s]
      if (h <= ncol(Bhat)) lab_cols[h] <- colnames(data$R)[s]
    }
    colnames(Bhat) <- ifelse(is.na(lab_cols),
                             paste0("spurious_", seq_along(lab_cols)), lab_cols)
    fp_Betas <- compute_FP_Betas(data$Betas, Bhat, tol = beta_tol)
    sign_Betas <- compute_sign_Betas(data$Betas, Bhat)

    Lam_in <- lambda_on_bins(est, bsp$train_bins)
    Lam_out <- lambda_on_bins(est, bsp$test_bins)
    rmse_lambda_in <- sqrt(mean(rowMeans(Lam_in - Lambda_true[bsp$train_bins, , drop = FALSE])^2))
    rmse_lambda_out <- sqrt(mean(rowMeans(Lam_out - Lambda_true[bsp$test_bins, , drop = FALSE])^2))
    rmse_counts_in <- sqrt(mean((Lam_in - mat_counts[bsp$train_bins, , drop = FALSE])^2))
    rmse_counts_out <- sqrt(mean((Lam_out - mat_counts[bsp$test_bins, , drop = FALSE])^2))

    data.frame(model = model_name,
               Kest = ncol(est$R_hat),
               Sensitivity = unname(sp["Sensitivity"]),
               Precision = unname(sp["Precision"]),
               F1 = unname(sp["F1"]),
               cosine_R = cosine_R,
               rmse_sig = rmse_sig,
               rmse_theta = rmse_theta,
               rmse_Betas = rmse_Betas,
               sign_Betas = sign_Betas,
               fp_Betas = fp_Betas,
               rmse_lambda_in = rmse_lambda_in,
               rmse_lambda_out = rmse_lambda_out,
               rmse_counts_in = rmse_counts_in,
               rmse_counts_out = rmse_counts_out,
               time = fit_runtime_mins(res),
               t(sampling_details_of(res, keep = colnames(est$R_hat))),
               stringsAsFactors = FALSE)
  }

  results <- rbind(
    eval_one("output_CompNMF.rds.gzip", "CompNMF", "CompNMF"),
    eval_one("output_SignatureAnalyzer.rds.gzip", "SignatureAnalyzer", "SignatureAnalyzer"),
    eval_one("output_map_FullModel.rds.gzip", "PPF", "PPF_map"),
    eval_one("output_mcmc_FullModel.rds.gzip", "PPF", "PPF_mcmc")
  )
  if (!is.null(results) && nrow(results) > 0) {
    results$Simulation <- basename(out_dir)
    results$Scenario <- basename(dirname(out_dir))
    results$n <- length(data$gr_Mutations)
    results$J <- data$simulation_parameters$J
  }
  results
}




#' The baseline of a fit, whichever package produced it
#'
#' SignaturePPF fits the activity prior, so the baseline is `$Baseline` and
#' `$Thetas` is the total. CompNMF has no covariates and stores the total in
#' `$Theta`; the caller divides that by the exposure itself.
ppf_baseline <- function(fit) {
  if (!is.null(fit$Baseline)) fit$Baseline else fit$Thetas
}


#' FIXED-signature evaluation of one fit: counts RMSE and probability calibration
postProcessCalibration_misspec <- function(out_dir,
                                           fit_file = "output_map_Fixed.rds.gzip",
                                           nbins = 10, beta_tol = 0.01,
                                           n_test_bins = 4000) {
  data <- readRDS(file.path(out_dir, "data.rds.gzip"))
  fit <- open_rds_file(file.path(out_dir, fit_file))
  if (is.null(fit) || is.null(data$p_true)) return(NULL)
  tilewidth <- data$simulation_parameters$tilewidth
  SignalTrackSim <- as.matrix(GenomicRanges::mcols(data$gr_SignalTrack))[, -1]
  CopyTrack_obs <- tilewidth * as.matrix(GenomicRanges::mcols(data$gr_CopyTrack)) / 2
  true_probs <- data$p_true
  true_sigs <- colnames(true_probs)

  bsp <- split_bins_misspec(data, n_test_bins = n_test_bins)

  n_bins <- nrow(SignalTrackSim)
  mat_counts <- unclass(table(factor(bsp$bin_mut, levels = seq_len(n_bins)),
                              factor(as.character(data$gr_Mutations$sample),
                                     levels = colnames(CopyTrack_obs))))
  storage.mode(mat_counts) <- "double"

  if (!is.null(fit$Betas)) {
    # ---- SignaturePPF (covariates + copy number) ----
    R <- fit$Signatures
    Theta_bl <- ppf_baseline(fit)
    Betas <- fit$Betas
    Lambda_hat <- reconstruct_lambda_raw(
      SignalTrack = SignalTrackSim, CopyTrack = CopyTrack_obs,
      Phi = Theta_bl, Betas = Betas)
    pred_probs <- Compute_mutation_Probs(data$gr_Mutations, R, Theta_bl, Betas)
  } else {
    # ---- CompNMF. The per-copy baseline is normalised over the TRAINING tiles
    #      only, so the held-out tiles are extrapolated by copy number alone.
    R <- fit$Signatures; Theta_tot <- fit$Theta
    Znorm <- colSums(CopyTrack_obs[bsp$train_bins, , drop = FALSE])
    denom <- t(Znorm)[rep(1, nrow(Theta_tot)), , drop = FALSE]
    Theta_bl <- Theta_tot / denom
    Beta0 <- matrix(0, ncol(SignalTrackSim), nrow(Theta_bl),
                    dimnames = list(colnames(SignalTrackSim), rownames(Theta_bl)))
    Lambda_hat <- reconstruct_lambda_raw(
      SignalTrack = SignalTrackSim, CopyTrack = CopyTrack_obs,
      Phi = Theta_bl, Betas = Beta0)
    pred_probs <- Compute_mutation_Probs(data$gr_Mutations, R, Theta_tot, Betas = NULL)
  }

  rmse_counts_in <- sqrt(mean((Lambda_hat[bsp$train_bins, , drop = FALSE] -
                                 mat_counts[bsp$train_bins, , drop = FALSE])^2))
  rmse_counts_out <- sqrt(mean((Lambda_hat[bsp$test_bins, , drop = FALSE] -
                                  mat_counts[bsp$test_bins, , drop = FALSE])^2))

  # The signature columns already carry catalogue names, so the true B is padded
  # onto the fitted shape - zeros for the distractors, which truly have no effect.
  Beta_hat <- if (!is.null(fit$Betas)) fit$Betas else
    matrix(0, ncol(SignalTrackSim), ncol(R),
           dimnames = list(colnames(SignalTrackSim), colnames(R)))
  Beta_true_pad <- pad_to_match(data$Betas, Beta_hat)
  rmse_Betas <- sqrt(mean((Beta_hat - Beta_true_pad)^2))
  sel <- beta_CI_excludes_zero(fit, dimnames_ref = dimnames(Beta_hat))
  fp_Betas <- compute_FP_Betas(Beta_true_pad, Beta_hat, tol = 0.02, selected = sel)
  sign_Betas <- compute_sign_Betas(Beta_true_pad, Beta_hat, selected = sel)

  lab <- colnames(pred_probs); lab[!(lab %in% true_sigs)] <- NA
  cal_on <- function(idx) evaluate_attribution_calibration(
    pred_probs[idx, , drop = FALSE], lab, data$signature_true[idx], true_sigs,
    nbins = nbins, p_oracle = true_probs[idx, , drop = FALSE])

  ac_in <- cal_on(bsp$train_mut)
  ac_out <- cal_on(bsp$test_mut)
  data.frame(model = sub("output_(.*)\\.rds\\.gzip$", "\\1", fit_file),
             rmse_counts_in = rmse_counts_in,
             rmse_counts_out = rmse_counts_out,
             rmse_Betas = rmse_Betas,
             sign_Betas = sign_Betas,
             fp_Betas = fp_Betas,
             attribution_acc_in = ac_in$accuracy,
             attribution_acc_out = ac_out$accuracy,
             estim_error_in = ac_in$estim_error,
             estim_error_out = ac_out$estim_error,
             brier_in = ac_in$brier, brier_out = ac_out$brier,
             ece_in = ac_in$ece, ece_out = ac_out$ece,
             mce_in = ac_in$mce, mce_out = ac_out$mce,
             oracle_l2_in = ac_in$oracle_l2, oracle_l2_out = ac_out$oracle_l2,
             oracle_kl_in = ac_in$oracle_kl, oracle_kl_out = ac_out$oracle_kl,
             n_mut_in = length(bsp$train_mut), n_mut_out = length(bsp$test_mut),
             Simulation = basename(out_dir), Scenario = basename(dirname(out_dir)),
             stringsAsFactors = FALSE)
}


#' Calibration curves: cumulative confidence against cumulative P(correct)
CALIBRATION_TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
CALIBRATION_COLS <- c("C>A" = "#16BDEB", "C>G" = "#000000", "C>T" = "#E22926",
                      "T>A" = "#A6A6A6", "T>C" = "#A1CE63", "T>G" = "#EBC6C4")


#' The cumulative calibration curve of one fit, split by substitution type
#'
#' Split out of `plot_calibration_curves` so that several fits can be assembled
#' into one faceted figure without duplicating the curve logic. `thin` keeps at
#' most that many equally spaced points per type, always including the last one:
#' the curve is a cumulative average over ~2e5 mutations, so drawing every point
#' is invisible on the page and expensive in the PDF.
calibration_curve_data <- function(data, fit, set = c("all", "in", "out"),
                                   n_test_bins = 4000, thin = NULL) {
  set <- match.arg(set)
  true_probs <- data$p_true
  pred_probs <- if (!is.null(fit$Betas))
    Compute_mutation_Probs(data$gr_Mutations, fit$Signatures, ppf_baseline(fit), fit$Betas)
  else
    Compute_mutation_Probs(data$gr_Mutations, fit$Signatures, fit$Theta, Betas = NULL)
  if (set != "all") {
    bsp <- split_bins_misspec(data, n_test_bins = n_test_bins)
    idx <- if (set == "in") bsp$train_mut else bsp$test_mut
    pred_probs <- pred_probs[idx, , drop = FALSE]
    true_probs <- true_probs[idx, , drop = FALSE]
  }
  n <- nrow(pred_probs)

  pred_idx <- max.col(pred_probs, ties.method = "first")
  best_prob <- pred_probs[cbind(seq_len(n), pred_idx)]
  pred_label <- colnames(pred_probs)[pred_idx]
  # True probability of the predicted signature: zero when it is a distractor
  # absent from the truth, i.e. a confidently wrong call.
  col_in_true <- match(pred_label, colnames(true_probs))
  p_correct <- numeric(n)
  ok <- !is.na(col_in_true)
  p_correct[ok] <- true_probs[cbind(which(ok), col_in_true[ok])]

  mut_type <- sub(".*\\[(.*)\\].*", "\\1", rownames(pred_probs))

  df <- do.call(rbind, lapply(CALIBRATION_TYPES, function(ty) {
    idx <- which(mut_type == ty)
    if (!length(idx)) return(NULL)
    m <- length(idx)
    conf <- cumsum(best_prob[idx]) / m
    pcorrect <- cumsum(p_correct[idx]) / m
    keep <- if (is.null(thin) || m <= thin) seq_len(m) else
      unique(c(round(seq(1, m, length.out = thin)), m))
    data.frame(type = ty, conf = conf[keep], pcorrect = pcorrect[keep],
               stringsAsFactors = FALSE)
  }))
  df$type <- factor(df$type, levels = CALIBRATION_TYPES)
  df
}


#' The endpoint of every curve, i.e. the overall mean of each substitution type
calibration_curve_ends <- function(df, by = "type") {
  do.call(rbind, lapply(split(df, df[by], drop = TRUE),
                        function(d) d[nrow(d), , drop = FALSE]))
}


plot_calibration_curves <- function(data, fit, set = c("all", "in", "out"),
                                    n_test_bins = 4000) {
  df <- calibration_curve_data(data, fit, set = set, n_test_bins = n_test_bins)
  ends <- calibration_curve_ends(df)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$conf, y = .data$pcorrect,
                                   colour = .data$type)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(data = ends, size = 2.4) +
    ggplot2::scale_colour_manual(values = CALIBRATION_COLS, name = "Mutation type") +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "Cumulative mean confidence", y = "Cumulative P(correct)") +
    ggplot2::theme_bw()
}


