## Produces: no figure. Helper file, sourced by the analysis scripts.

## Goodness-of-fit and patient-level deviation diagnostics for the PPF model.
##
## Two referee points are answered here.
##
##   (i)  MODEL CHECKING. Does the fitted intensity reproduce the observed counts
##        at 2 kb, and are there systematic residuals by genomic region, by
##        chromosome, or by mutation type?
##   (ii) PATIENT-LEVEL HETEROGENEITY. Where do individual patients depart from
##        the model, and is that departure driven by mutation burden?
##
## Nothing here is a posterior predictive check. Under the model the counts in
## disjoint bins are INDEPENDENT POISSON with a known mean, so every quantity a
## predictive simulation would estimate is available in closed form - including
## the sampling variability the check has to be judged against. Simulating from
## the fit would only add Monte Carlo error to numbers we can write down.
##
## ---------------------------------------------------------------------------
## THE LARGE-n PROBLEM, AND WHAT IS DONE ABOUT IT
##
## The cohort has 707k mutations on 1.39M bins. Any consistent test - chi-square,
## Kolmogorov-Smirnov on the rescaled event times - rejects at that sample size
## for deviations far below the size that would change a conclusion. A p-value is
## therefore not a measurement here, and the diagnostics are built around three
## substitutes.
##
##   1. EFFECT SIZES, not p-values. The dispersion ratio D = X^2 / R estimates a
##      fixed quantity and does not grow with n. If the model's regional rate
##      carries a multiplicative error of coefficient of variation `cv` on top of
##      Poisson noise, then
##
##          E[X^2] = sum_r (e_r + cv^2 e_r^2) / e_r = R + cv^2 sum_r e_r,
##          so    D = 1 + cv^2 * mean(e),   i.e.   cv = sqrt((D - 1) / mean(e)).
##
##      `cv` is reported alongside D: it is the fractional error in the predicted
##      regional rate, beyond Poisson noise, and it is what a reader can act on.
##      "The 1 Mb rate is right to within 8%" is a statement about the model;
##      "p < 1e-300" is a statement about the sample size.
##
##      D is judged against its EXACT Poisson sampling band rather than against
##      chi-square asymptotics, so it stays valid at 2 kb where e_r << 1. For
##      X ~ Pois(e) and Y = (X - e)^2 / e, E[Y] = 1 and Var(Y) = 2 + 1/e, hence
##
##          sd(D) = sqrt(sum_r (2 + 1/e_r)) / R.
##
##      No minimum expected count is imposed anywhere: the identity above holds
##      for any e_r > 0, and filtering on e would select regions on a quantity
##      correlated with the residual being tested.
##
##   2. EQUAL POWER, where a test is still wanted. The time-rescaling test is run
##      after binomially thinning each slice to a common expected number of
##      events n0. Thinning a Poisson process with probability p leaves a Poisson
##      process of intensity p*lambda, so the test remains EXACT; what changes is
##      that a 64k-mutation hypermutator and a 1.2k-mutation patient are now
##      judged at the same power. Without it the per-patient p-values order the
##      patients by burden and say nothing about fit.
##
##   3. LOCALISATION instead of rejection. X^2 is decomposed across regions and
##      the concentration of the excess is reported: what share of it sits in the
##      top 0.1% of regions, and what D falls to once they are removed. A model
##      whose excess dispersion is carried by a handful of hotspots is not the
##      same object as one that is wrong everywhere, and the referee's question
##      ("systematic residuals in specific genomic regions") is exactly the
##      question of which of the two this is.
##
## ---------------------------------------------------------------------------
## WHAT THE INTENSITY IS, AND WHY THE CHANNEL SPLIT IS EXACT
##
## The model's expected count in bin b for patient j is
##
##     lambda_bj = c_j(b) * sum_k phi_kj * exp(beta_k' x_b),
##
## the sum over the 96 channels having already collapsed because the signature
## columns are normalised. For a SUBSET S of channels the same reconstruction
## holds with each signature weighted by the mass it puts on S:
##
##     lambda^S_bj = c_j(b) * sum_k w_k(S) * phi_kj * exp(beta_k' x_b),
##     w_k(S) = sum_{i in S} r_ik.
##
## This is exact, not an approximation - the intensity is linear in the signature
## matrix - so the six macro mutation classes are six reweightings of one
## precomputed exp(X beta), and cost one matrix-vector product each.
##
## Pooling over patients is equally exact: the superposition of independent
## Poisson processes is a Poisson process whose intensity is the sum. Every
## pooled slice below is that superposition.
##
## ---------------------------------------------------------------------------
## Usage
##
##   g <- gof_setup(data, fit)                        # one pass, then reuse
##   gof_count_distribution(g)                        # 2 kb marginal, closed form
##   gof_dispersion_scales(g)                         # D and cv, class x scale
##   reg <- gof_regions(g, width = 1e6)               # residuals along the genome
##   pat <- gof_patients(g)                           # per-patient deviations
##
## The fit must carry `Signatures`, `Baseline` and `Betas`. `Baseline` is phi,
## NOT `Thetas`: under the activity prior theta = phi * q_j(beta_k), and passing
## Thetas would inflate every intensity by q_j(beta_k).


## ------------------------------------------------------------------ constants

MUT_CLASSES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
MUT_COLS <- c("C>A" = "#16BDEB", "C>G" = "#000000", "C>T" = "#E22926",
              "T>A" = "#A6A6A6", "T>C" = "#A1CE63", "T>G" = "#EBC6C4",
              "all" = "#000D8B")

## The macro class of a trinucleotide channel: "A[C>T]G" -> "C>T".
channel_class <- function(x) sub(".*\\[(.*)\\].*", "\\1", as.character(x))


#' Region index of every bin, at a given region width
#'
#' Regions are cut on genomic coordinates and never straddle a chromosome, which
#' a running index over bin position (`(bin - 1) %/% n + 1`) would do at every
#' chromosome boundary. The bins arrive in genomic order with each chromosome
#' contiguous, so the region id is one cumulative sum.
#'
#' @param chrom,start Chromosome and start coordinate of every bin.
#' @param width Region width in base pairs.
#' @return An integer vector of region ids, consecutive from 1 in genomic order.
region_index <- function(chrom, start, width) {
  blk <- floor((as.numeric(start) - 1) / width)
  n <- length(blk)
  if (n == 0L) return(integer(0))
  brk <- c(TRUE, blk[-1] != blk[-n] | chrom[-1] != chrom[-n])
  cumsum(brk)
}


#' Everything the diagnostics need, computed once
#'
#' The expensive part is `exp(X beta)`, an `n_bins x K` matrix, and its product
#' with the copy-number track. Both are formed here and reused by every
#' diagnostic, because at 2 kb they cost more than all the tests together.
#'
#' `EA[b, k] = exp(beta_k' x_b) * sum_j c_j(b) phi_kj` is the per-bin, per-
#' signature intensity ALREADY SUMMED OVER PATIENTS, so a pooled intensity for
#' any subset of channels is a single matrix-vector product `EA %*% w`. The
#' per-patient intensities need `E` and `Phi` separately and are formed on demand,
#' one patient at a time, rather than materialising an `n_bins x n_samples`
#' matrix that would be 1.2 GB at 2 kb.
#'
#' @param data A cohort object from [load_cohort()], carrying `gr_SignalTrack`.
#' @param fit A `SignaturePPF` fit with `Signatures`, `Baseline` and `Betas`.
#' @param bin_of_mut Bin index of each mutation. Computed if not supplied.
#' @param widths Region widths, in base pairs, to precompute region indices for.
#' @return A list consumed by every other function in this file.
gof_setup <- function(data, fit, bin_of_mut = NULL,
                      widths = c(1e4, 1e5, 1e6, 1e7)) {
  stopifnot(all(c("Signatures", "Baseline", "Betas") %in% names(fit)))
  if (is.null(data$gr_SignalTrack))
    stop("`data` has no `gr_SignalTrack`; load it with load_cohort().", call. = FALSE)

  samples <- intersect(colnames(fit$Baseline), colnames(data$CopyTrack))
  if (!length(samples)) stop("the fit and the cohort share no samples.", call. = FALSE)
  Phi <- fit$Baseline[, samples, drop = FALSE]
  Copy <- data$CopyTrack[, samples, drop = FALSE]

  ## The same clamp reconstruct_lambda() uses, so the diagnostics score exactly
  ## the intensity the rest of the pipeline reports.
  E <- exp(pmin(pmax(data$SignalTrack[, rownames(fit$Betas), drop = FALSE] %*%
                       fit$Betas, -20), 20))
  ## sum_j c_j(b) phi_kj, then folded into E: the patient-pooled intensity of
  ## signature k in bin b, per unit of signature mass.
  EA <- E * (Copy %*% t(Phi))

  ## Class weights: w[k, s] is the mass signature k puts on macro class s.
  cls_of_channel <- channel_class(rownames(fit$Signatures))
  W <- vapply(MUT_CLASSES, function(s)
    colSums(fit$Signatures[cls_of_channel == s, , drop = FALSE]),
    numeric(ncol(fit$Signatures)))
  W <- cbind(W, all = rowSums(W))

  ## Mutations: keep only those inside a bin and belonging to a fitted sample.
  if (is.null(bin_of_mut)) bin_of_mut <- bin_of_mutation(data)
  smp <- as.character(data$gr_Mutations$sample)
  keep <- !is.na(bin_of_mut) & smp %in% samples
  if (any(!keep))
    message(sum(!keep), " of ", length(keep),
            " mutations are outside every bin or belong to an unfitted sample.")

  gr <- data$gr_SignalTrack
  chrom <- as.character(GenomicRanges::seqnames(gr))
  regions <- lapply(widths, function(w)
    region_index(chrom, GenomicRanges::start(gr), w))
  names(regions) <- format(widths, scientific = FALSE, trim = TRUE)

  list(E = E, Phi = Phi, EA = EA, CopyTrack = Copy, W = W,
       n_bins = nrow(E), samples = samples,
       bin_of_mut = bin_of_mut[keep],
       sample_of_mut = smp[keep],
       class_of_mut = channel_class(data$gr_Mutations$channel)[keep],
       chrom = chrom,
       start = GenomicRanges::start(gr),
       end = GenomicRanges::end(gr),
       bin_weight = if (!is.null(gr$bin_weight)) as.numeric(gr$bin_weight)
                    else GenomicRanges::width(gr),
       regions = regions, widths = widths)
}


#' Expected counts per bin under the fit
#'
#' @param g Output of [gof_setup()].
#' @param class One of `MUT_CLASSES`, or `"all"`.
#' @param sample A single sample id, or `NULL` for the patient-pooled intensity.
#' @return A numeric vector of length `g$n_bins`.
gof_lambda <- function(g, class = "all", sample = NULL) {
  w <- g$W[, class]
  if (is.null(sample)) return(as.numeric(g$EA %*% w))
  as.numeric(g$CopyTrack[, sample] * (g$E %*% (w * g$Phi[, sample])))
}


#' Observed counts per bin
#'
#' @inheritParams gof_lambda
#' @return An integer vector of length `g$n_bins`.
gof_obs <- function(g, class = "all", sample = NULL) {
  sel <- rep(TRUE, length(g$bin_of_mut))
  if (class != "all") sel <- sel & g$class_of_mut == class
  if (!is.null(sample)) sel <- sel & g$sample_of_mut == sample
  tabulate(g$bin_of_mut[sel], nbins = g$n_bins)
}


## ---------------------------------------------------------------------------
## 1. The 2 kb marginal distribution of counts, in closed form.
##
##    The referee asks whether the model captures the observed distribution of
##    counts at 2 kb. Under the model the count in cell (bin, patient) is
##    Poisson(lambda_bj) independently, so the expected NUMBER OF CELLS holding
##    exactly m mutations is
##
##        N_m = sum_{b,j} Pois(m; lambda_bj),
##
##    and because the cells are independent Bernoulli indicators, its exact
##    variance is sum_{b,j} p_m (1 - p_m). Observed against expected, with a
##    z-score that needs no simulation and no asymptotics.
##
##    The interpretable column is `ratio` = observed / expected. At 2 kb the mean
##    count is ~0.005, so the model predicts very few multi-mutation bins; an
##    excess there is exactly the local clustering (kataegis, hotspots) that a
##    smooth intensity cannot express, and the ratio says how large it is.
## ---------------------------------------------------------------------------

#' Observed and model-implied distribution of per-bin counts
#'
#' `ratio` is the headline: the factor by which the model under- or over-predicts
#' the number of bins holding exactly that many mutations. `mut_share` says how
#' much of the cohort those bins actually hold, which is what keeps a large ratio
#' in proportion - the model can be out by an order of magnitude on the bins
#' carrying five mutations and still be describing 99% of the data.
#'
#' @param g Output of [gof_setup()].
#' @param class One of `MUT_CLASSES`, or `"all"`.
#' @param mmax Counts are tabulated for `0, 1, ..., mmax` and `>= mmax + 1`.
#' @return A data frame, one row per count value.
gof_count_distribution <- function(g, class = "all", mmax = 4L) {
  ms <- seq.int(0L, mmax)

  ## ---- expected, accumulated one patient at a time
  exp_n <- numeric(length(ms) + 1L)   # counts 0..mmax, then the tail
  var_n <- numeric(length(ms) + 1L)
  exp_total <- 0                      # expected mutations, for the tail row
  w <- g$W[, class]
  for (j in g$samples) {
    lam <- as.numeric(g$CopyTrack[, j] * (g$E %*% (w * g$Phi[, j])))
    exp_total <- exp_total + sum(lam)
    p <- exp(-lam)                    # the one exp() per cell; the rest recurses
    tail <- rep(1, length(lam))
    for (i in seq_along(ms)) {
      if (ms[i] > 0L) p <- p * lam / ms[i]
      exp_n[i] <- exp_n[i] + sum(p)
      var_n[i] <- var_n[i] + sum(p * (1 - p))
      tail <- tail - p
    }
    tail <- pmax(tail, 0)             # rounding only; the true value is >= 0
    exp_n[length(exp_n)] <- exp_n[length(exp_n)] + sum(tail)
    var_n[length(var_n)] <- var_n[length(var_n)] + sum(tail * (1 - tail))
  }

  ## ---- observed. Only cells with at least one mutation need to be touched:
  ##      everything else is the zero cell, by subtraction.
  sel <- if (class == "all") rep(TRUE, length(g$bin_of_mut))
         else g$class_of_mut == class
  cell <- table(paste(g$bin_of_mut[sel], g$sample_of_mut[sel]))
  occ <- tabulate(as.integer(cell), nbins = mmax + 1L)
  n_cells <- as.numeric(g$n_bins) * length(g$samples)
  obs_n <- c(n_cells - length(cell), occ[seq_len(mmax)],
             length(cell) - sum(occ[seq_len(mmax)]))

  ## Mutations held by each group of cells. Exact for the counted values; the
  ## tail follows by subtraction from the totals, observed and expected alike.
  n_mut <- sum(sel)
  obs_mut <- c(0, ms[-1] * obs_n[seq(2, length(ms))],
               n_mut - sum(ms * obs_n[seq_along(ms)]))
  exp_mut <- c(0, ms[-1] * exp_n[seq(2, length(ms))],
               exp_total - sum(ms * exp_n[seq_along(ms)]))

  data.frame(class = class,
             count = c(as.character(ms), paste0(">=", mmax + 1L)),
             observed = obs_n, expected = exp_n,
             ratio = obs_n / exp_n,
             z = (obs_n - exp_n) / sqrt(pmax(var_n, .Machine$double.eps)),
             mutations = obs_mut, mut_share = obs_mut / n_mut,
             mutations_expected = exp_mut,
             row.names = NULL, stringsAsFactors = FALSE)
}


## ---------------------------------------------------------------------------
## 2. Dispersion as an effect size, across aggregation scales.
## ---------------------------------------------------------------------------

#' Aggregate a per-bin vector onto regions
#'
#' `rowsum(reorder = FALSE)` keeps first-appearance order, which for the ids
#' [region_index()] produces is genomic order.
agg_regions <- function(x, region) as.numeric(rowsum(x, region, reorder = FALSE))


#' Dispersion of the observed counts about the fitted intensity
#'
#' @param obs,exp Observed and expected counts, already aggregated onto regions.
#' @param trim Fraction of regions, ordered by their contribution to `X^2`, to
#'   drop for the trimmed statistics. The default 0.001 is one region in a
#'   thousand.
#' @return A one-row data frame; see the file header for what `D`, `cv_excess`
#'   and `sd_D` mean and why they are preferred to a p-value.
gof_dispersion <- function(obs, exp, trim = 0.001) {
  keep <- is.finite(exp) & exp > 0
  o <- obs[keep]; e <- exp[keep]
  R <- length(e)
  if (R < 2L) return(NULL)

  contrib <- (o - e)^2 / e
  X2 <- sum(contrib)
  D <- X2 / R
  ## Exact Poisson sampling sd of D: Var((X - e)^2 / e) = 2 + 1/e.
  sd_D <- sqrt(sum(2 + 1 / e)) / R

  ord <- order(contrib, decreasing = TRUE)
  n_top <- max(1L, round(trim * R))
  top <- ord[seq_len(n_top)]
  X2_trim <- X2 - sum(contrib[top])

  data.frame(
    n_regions = R,
    obs_total = sum(o), exp_total = sum(e),
    mean_exp = mean(e),
    D = D, sd_D = sd_D,
    z_D = (D - 1) / sd_D,
    ## The fractional error in the regional rate implied by D, beyond Poisson.
    cv_excess = sqrt(pmax(D - 1, 0) / mean(e)),
    trim = trim,
    D_trim = X2_trim / (R - n_top),
    cv_excess_trim = sqrt(pmax(X2_trim / (R - n_top) - 1, 0) / mean(e[-top])),
    ## How concentrated the excess is: the share of the total X^2 carried by the
    ## trimmed regions, and by the single worst one.
    top_share = sum(contrib[top]) / X2,
    max_share = max(contrib) / X2,
    row.names = NULL, stringsAsFactors = FALSE)
}


## ---------------------------------------------------------------------------
## 3. Residuals along the genome.
## ---------------------------------------------------------------------------

#' Exact two-sided Poisson p-value, R's convention
#'
#' The p-value `stats::poisson.test(x, T = 1, r = e)` returns: the total
#' probability of every count no more likely than the one observed. That is NOT
#' the same as "twice the smaller tail" on a discrete asymmetric distribution -
#' on this cohort the two disagree on 21 of 2,867 megabases - and this is the
#' convention a referee can reproduce with one line, so it is the one used.
#'
#' `poisson.test()` is not vectorised and builds an htest object per call, which
#' is why it is wrapped rather than called directly at the (patient, region)
#' scale; the wrapper is still a loop, costing about 30 s for 324k cells, but it
#' is provably the same number.
#'
#' WHY NOT CHI-SQUARE. At the pooled megabase scale it would do: the expected
#' counts are ~240 and the two tests disagree on 0.7% of regions. At the
#' (patient, megabase) scale the median expected count is 1.2 and 97% of cells
#' are below 10, where the chi-square approximation is badly anti-conservative -
#' in a quarter of cells it returns less than half the exact p-value. One test
#' that is correct at both scales is cheaper than two conventions to explain.
#'
#' @param obs,exp Observed counts and fitted means, of equal length.
#' @return A numeric vector of p-values.
poisson_p_exact <- function(obs, exp) {
  mapply(function(x, r) {
    if (!is.finite(r) || r <= 0) return(NA_real_)
    stats::poisson.test(as.integer(x), T = 1, r = r)$p.value
  }, obs, exp)
}


#' Observed, expected and Pearson residual per genomic region
#'
#' The p-value is the EXACT two-sided Poisson tail probability, not the normal
#' approximation to the Pearson residual: for a single patient or a single
#' mutation class the expected regional count is small, and the approximation is
#' unreliable exactly there.
#'
#' A flagged region here is a claim about a REGION, not about the model as a
#' whole - which is why the multiplicity adjustment is over the regions tested
#' and the interesting summary is how few of them there are, and how much of the
#' total Pearson chi-square they carry.
#'
#' @param g Output of [gof_setup()].
#' @param width Region width in base pairs; must be one of `g$widths`.
#' @param class One of `MUT_CLASSES`, or `"all"`.
#' @param sample A single sample id, or `NULL` to pool over patients.
#' @param p_adjust Multiplicity adjustment, passed to [stats::p.adjust()].
#' @return A data frame, one row per region, in genomic order.
gof_regions <- function(g, width = 1e6, class = "all", sample = NULL,
                        p_adjust = "BH") {
  key <- format(width, scientific = FALSE, trim = TRUE)
  reg <- g$regions[[key]]
  if (is.null(reg))
    stop("no region index for width ", key,
         ". Rebuild with gof_setup(widths = ...).", call. = FALSE)

  lam <- agg_regions(gof_lambda(g, class, sample), reg)
  obs <- agg_regions(gof_obs(g, class, sample), reg)
  first <- !duplicated(reg)
  last <- !duplicated(reg, fromLast = TRUE)

  pv <- poisson_p_exact(obs, lam)
  contrib <- (obs - lam)^2 / pmax(lam, .Machine$double.eps)

  ## Usable sequence actually modelled in the window, as a fraction of its span.
  ## The exposure is already proportional to it, so it cannot bias the residual -
  ## it is here so a reader can confirm that, rather than having to trust it.
  usable <- agg_regions(g$bin_weight, reg)
  span <- agg_regions(as.numeric(g$end - g$start + 1), reg)

  out <- data.frame(
    region = seq_along(obs),
    chrom = g$chrom[first], start = g$start[first], end = g$end[last],
    n_bins = as.integer(agg_regions(rep(1L, g$n_bins), reg)),
    usable_frac = usable / span,
    observed = obs, expected = lam,
    resid = (obs - lam) / sqrt(pmax(lam, .Machine$double.eps)),
    contrib = contrib, contrib_share = contrib / sum(contrib),
    p_val = pv, row.names = NULL, stringsAsFactors = FALSE)
  out$padj <- stats::p.adjust(out$p_val, method = p_adjust)
  out$class <- class
  out
}


#' Locus class of a set of windows
#'
#' Only the immunoglobulin loci. Assembly gaps and the ENCODE blacklist were
#' TRIED here and removed again, because on this cohort they annotate nothing:
#' the preprocessing already subtracts masked sequence from each bin's usable
#' width, and zero of the 707,104 mutations fall inside either mask. A megabase
#' "overlapping the blacklist" merely touches a masked interval somewhere in its
#' million bases while being 99%+ callable, so the overlap carries no information
#' about why its residual is large. Use `usable_frac` from [gof_regions()] if you
#' want to check a window for masking - that is the quantity that would matter,
#' and it is ~100% everywhere the residuals are largest.
#'
#' The immunoglobulin loci are different: their mutations ARE in the data, and
#' somatic hypermutation there is biology the intensity model never claimed to
#' describe, so an excess is expected rather than a failure.
#'
#' @param reg A data frame with `chrom`, `start`, `end`.
#' @return A character vector: the locus, or `""`.
annotate_regions_class <- function(reg) {
  if (!nrow(reg)) return(character(0))
  gr <- GenomicRanges::GRanges(reg$chrom,
                               IRanges::IRanges(reg$start, reg$end))
  ## Canonical hg19 extents. Rounded-outward bounds are not harmless here: a
  ## 1.7 Mb "IGH" reaching to 107.3 Mb swallows the subtelomeric tip of chr14,
  ## and labelled the single worst window in the genome as immunoglobulin when it
  ## lies 120 kb past the end of the locus.
  ig <- GenomicRanges::GRanges(
    c("chr14", "chr2", "chr22"),
    IRanges::IRanges(c(105586437, 88857361, 22380474),
                     c(106879844, 90235368, 23265085)))
  ig$locus <- c("IGH", "IGK", "IGL")
  out <- rep("", nrow(reg))
  h <- suppressWarnings(GenomicRanges::findOverlaps(gr, ig))
  if (length(h))
    out[S4Vectors::queryHits(h)] <- ig$locus[S4Vectors::subjectHits(h)]
  out
}


#' Coordinate labels
#'
#' `chr14:106 Mb`, and nothing else. Every richer label tried here has been a
#' claim the data did not support - gene symbols read as "the model fails at this
#' gene" when they are only the longest gene in the window, and a locus name
#' asserts a mechanism that has to hold for that window rather than for a range
#' it happens to sit in. The locus is still recorded in the tables, where a
#' reader can weigh it; the figure states the coordinate.
region_labels <- function(reg) {
  sprintf("%s:%.0f Mb", reg$chrom, reg$start / 1e6)
}


## ---------------------------------------------------------------------------
## 4. Time-rescaling test, at equal power.
##
##    Under the random time change theorem, mapping the events of a Poisson
##    process through its cumulative intensity gives a rate-1 homogeneous
##    process: conditional on their number the rescaled positions are iid
##    U(0, 1). That conditional form is the right one here, because the total is
##    a FITTED quantity - the model matches it by construction, so a test that
##    also charges for the total would spend its power on something already
##    known to agree. What is left is the regional SHAPE, which is what the
##    covariates are claimed to explain.
##
##    The fitted intensity is piecewise constant on a bin, so an event's position
##    inside its bin is uniform under the model and is drawn at random rather
##    than taken from the base pair (a randomised PIT; Dunn & Smyth 1996). This
##    removes the ties a bin-level position necessarily creates and points the
##    test at the BETWEEN-bin structure. Using the true base pair would also test
##    within-bin uniformity, which the model never claimed.
## ---------------------------------------------------------------------------

#' Kolmogorov-Smirnov distance from U(0, 1)
#'
#' Reported as an effect size at the full event count, where the p-value is
#' uninformative. Under a correct model it decays like 0.87 / sqrt(n), which is
#' the reference `ks_null` returns.
ks_distance <- function(u) {
  n <- length(u)
  if (n < 2L) return(NA_real_)
  u <- sort(u)
  max(max(seq_len(n) / n - u), max(u - (seq_len(n) - 1) / n))
}

#' Expected KS distance under the null, for reference against [ks_distance()]
ks_null <- function(n) sqrt(pi / 2) * log(2) / sqrt(n)

#' Asymptotic two-sided Kolmogorov-Smirnov p-value
#'
#' `stats::ks.test()` switches to exactly this series above n = 100, and every
#' patient here has at least a thousand events, so calling it directly is the
#' same number without the overhead of building a htest object per replicate.
ks_pvalue <- function(D, n) {
  t <- sqrt(n) * D
  if (!is.finite(t) || t <= 0) return(1)
  k <- 1:200
  min(1, max(0, 2 * sum((-1)^(k - 1) * exp(-2 * k^2 * t^2))))
}


#' Time-rescaling test on one slice, at full power and at a fixed power
#'
#' TWO tests are returned, because two different questions are asked of them and
#' one set of numbers cannot answer both.
#'
#'   `p_full`  the test on every event the patient has. This is the one that
#'             says whether a patient is misspecified, and it should be used
#'             for that: a patient with 64,000 mutations genuinely does carry
#'             more evidence than one with 1,200, and discarding it to make the
#'             two comparable discards real power.
#'   `p_thin`  the test after binomially thinning to `n0` events. Only for
#'             asking whether misfit is ASSOCIATED WITH BURDEN, where the
#'             burden-driven power difference is the confounder and has to go.
#'
#' Thinning does not shrink the deviation, it raises the bar: a spike of k events
#' moves the empirical CDF by about k/n, which thinning leaves alone, while the
#' critical distance grows like 1/sqrt(n). So a focal excess that is decisive at
#' full n can sit under the threshold after thinning. That is the correct
#' behaviour for a power-matched comparison and the wrong behaviour for
#' detection, hence both columns.
#'
#' The within-bin offset is redrawn on every replicate along with the thinning
#' (a randomised PIT; Dunn & Smyth 1996), and the reported values are medians
#' over replicates, so neither statistic rests on one random tie-break.
#'
#' @param lambda Expected counts per bin, in genomic order.
#' @param ev_bins Bin index of each event, into `lambda`.
#' @param n0 Events to retain for the thinned arm. `Inf` skips it.
#' @param n_rep Replicates.
#' @param min_events Slices with fewer events return `NA`.
gof_rescale_test <- function(lambda, ev_bins, n0 = Inf, n_rep = 25,
                             min_events = 50) {
  n <- length(ev_bins)
  out <- list(n = n, n_thin = NA_integer_, ks_full = NA_real_,
              p_full = NA_real_, ks_gap = NA_real_, p_gap = NA_real_,
              ks_null = NA_real_, p_thin = NA_real_,
              u = numeric(0), u_full = numeric(0), lam_total = NA_real_)
  tot <- sum(lambda)
  if (n < min_events || !is.finite(tot) || tot <= 0) return(out)

  Lam_left <- cumsum(lambda) - lambda
  keep_p <- if (is.finite(n0)) min(1, n0 / n) else 1

  ksf <- pf <- pt <- ksg <- pg <- numeric(n_rep); nk <- integer(n_rep)
  u_keep <- NULL; u_all_keep <- NULL
  for (r in seq_len(n_rep)) {
    u_all <- (Lam_left[ev_bins] + lambda[ev_bins] * runif(n)) / tot
    u_all <- pmin(pmax(u_all, 0), 1)
    ksf[r] <- ks_distance(u_all)
    pf[r] <- ks_pvalue(ksf[r], n)
    ## Brown et al. (2002): the rescaled GAPS are iid Exp(1), so
    ## z = 1 - exp(-Delta) is iid Uniform(0, 1). This is the statistic the KS and
    ## exponential panels actually draw, and it is a different test from the one
    ## on the rescaled positions above - sharper against local clustering.
    z <- 1 - exp(-diff(c(0, sort(u_all))) * tot)
    ksg[r] <- ks_distance(z)
    pg[r] <- ks_pvalue(ksg[r], n)

    u <- if (keep_p >= 1) u_all else u_all[runif(n) < keep_p]
    nk[r] <- length(u)
    pt[r] <- if (length(u) < 2L) NA_real_ else ks_pvalue(ks_distance(u), length(u))
    if (is.null(u_keep)) { u_keep <- sort(u); u_all_keep <- sort(u_all) }
  }
  list(n = n, n_thin = as.integer(round(mean(nk))),
       ks_full = stats::median(ksf), p_full = stats::median(pf),
       ks_gap = stats::median(ksg), p_gap = stats::median(pg),
       ks_null = ks_null(n), p_thin = stats::median(pt, na.rm = TRUE),
       u = if (is.null(u_keep)) numeric(0) else u_keep,
       u_full = if (is.null(u_all_keep)) numeric(0) else u_all_keep,
       lam_total = tot)
}


## ---------------------------------------------------------------------------
## 5. Patient-level deviations.
##
##    Not a hunt for patients to exclude. The per-patient total is fitted, so a
##    patient's overall residual is ~0 by construction and says nothing; what is
##    testable is whether the model gets the WITHIN-patient regional distribution
##    right, and whether it gets it less right for the hypermutated ones.
##
##    Three quantities per patient, each answering a different form of the
##    question:
##      D          how far the regional counts scatter beyond Poisson (effect size)
##      p_unif     whether the regional shape is rejected, at power held equal
##                 across patients by thinning
##      focal      what share of the patient's own misfit sits in the worst 0.1%
##                 of its regions - the difference between a patient the model
##                 describes badly everywhere and one with a few hotspots
## ---------------------------------------------------------------------------

#' Per-patient goodness of fit
#'
#' @param g Output of [gof_setup()].
#' @param width Region width for the dispersion statistics.
#' @param n0 Events retained per patient. `NULL` (default) sets it to the
#'   smallest burden in the cohort, which is the only value at which the test is
#'   ACTUALLY at equal power: a patient with fewer than `n0` events cannot be
#'   thinned up to it, so it keeps all of its events and is tested at lower power
#'   than the rest, while appearing in the same figure. An explicit `n0` above
#'   the minimum is allowed and warns, naming how many patients it affects.
#' @param n_rep Passed to [gof_rescale_test()].
#' @param seed Set once, so the randomised PIT and the thinning are reproducible.
#' @param qq Also return one equal-power QQ curve per patient.
#' @param n_points Points per QQ curve.
#' @param verbose Report progress; the loop is ~1 s a patient at 2 kb.
#' @return A data frame with one row per patient, or, if `qq = TRUE`, a list of
#'   that data frame and the curves.
gof_patients <- function(g, width = 1e6, n0 = NULL, n_rep = 25, seed = 1,
                         qq = FALSE, n_points = 200, verbose = TRUE) {
  set.seed(seed)
  key <- format(width, scientific = FALSE, trim = TRUE)
  reg <- g$regions[[key]]
  if (is.null(reg)) stop("no region index for width ", key, call. = FALSE)

  burden <- as.integer(table(factor(g$sample_of_mut, levels = g$samples)))
  if (is.null(n0)) {
    n0 <- min(burden)
    if (verbose) message("n0 = ", n0, " (the smallest burden in the cohort)")
  } else if (any(burden < n0)) {
    warning(sum(burden < n0), " of ", length(burden), " patients have fewer ",
            "than n0 = ", n0, " events (smallest ", min(burden), "), so they ",
            "are NOT thinned and are tested at lower power than the rest. ",
            "Set n0 <= ", min(burden), " for the comparison to be at equal ",
            "power.", call. = FALSE)
  }
  ## Genome actually modelled, for the burden-per-Mb column.
  mb <- sum(g$bin_weight) / 1e6

  rows <- vector("list", length(g$samples))
  curves <- vector("list", length(g$samples))
  for (i in seq_along(g$samples)) {
    j <- g$samples[i]
    if (verbose && i %% 20 == 0) message("  patient ", i, " / ", length(g$samples))
    lam <- gof_lambda(g, "all", j)
    ev <- g$bin_of_mut[g$sample_of_mut == j]

    d <- gof_dispersion(agg_regions(tabulate(ev, nbins = g$n_bins), reg),
                        agg_regions(lam, reg))
      ks <- gof_rescale_test(lam, ev, n0 = n0, n_rep = n_rep)

    rows[[i]] <- data.frame(
      sample = j, burden = length(ev), expected = sum(lam),
      burden_per_Mb = length(ev) / mb,
      D = d$D, sd_D = d$sd_D, z_D = d$z_D, cv_excess = d$cv_excess,
      D_trim = d$D_trim, focal = d$top_share,
      ks_full = ks$ks_full, ks_null = ks$ks_null,
      ks_ratio = ks$ks_full / ks$ks_null,
      ks_crit = ks_critical(length(ev)),
      p_full = ks$p_full,
      ks_gap = ks$ks_gap, p_gap = ks$p_gap,
      n_thin = ks$n_thin, p_thin = ks$p_thin,
      row.names = NULL, stringsAsFactors = FALSE)

    if (qq && length(ks$u) > 1L) {
      n <- length(ks$u)
      idx <- unique(round(seq(1, n, length.out = min(n_points, n))))
      curves[[i]] <- data.frame(sample = j, theoretical = (idx - 0.5) / n,
                                empirical = ks$u[idx],
                                row.names = NULL, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  ## Both arms are adjusted, and NEITHER is called just `padj`: which one a
  ## figure uses is a choice about what is being asked, and a bare name would
  ## let a call site make it by accident.
  out$padj_full <- stats::p.adjust(out$p_full, method = "BH")
  out$padj_gap <- stats::p.adjust(out$p_gap, method = "BH")
  out$padj_thin <- stats::p.adjust(out$p_thin, method = "BH")
  if (!qq) return(out)
  list(patients = out, qq = do.call(rbind, curves))
}


#' Observed and expected counts for every (region, patient) cell
#'
#' The unit the volcano is drawn in, and the unit that answers "is this region
#' flagged because the cohort departs from the model there, or because ONE
#' patient does". Both questions need the same table, so it is built once.
#'
#' 2,867 regions x 113 patients is 324k cells at 1 Mb, which is small; what is
#' expensive is the intensity, so the loop runs one patient at a time and
#' aggregates immediately rather than holding a bins-by-patients matrix.
#'
#' @param g Output of [gof_setup()].
#' @param width Region width in base pairs; must be one of `g$widths`.
#' @param class One of `MUT_CLASSES`, or `"all"`.
#' @param p_adjust Multiplicity adjustment over all cells.
#' @param compute_p Exact Poisson p-values are ~30 s for 324k cells; the
#'   patient-driver decomposition needs only the counts, so it skips them.
#' @return A data frame with one row per (region, patient).
gof_region_patient <- function(g, width = 1e6, class = "all", p_adjust = "BH",
                               compute_p = TRUE) {
  key <- format(width, scientific = FALSE, trim = TRUE)
  rid <- g$regions[[key]]
  if (is.null(rid)) stop("no region index for width ", key, call. = FALSE)
  R <- max(rid)

  ## Region of every mutation, once: the per-patient step is then a tabulate
  ## over R bins rather than over 1.39M.
  rom <- rid[g$bin_of_mut]
  first <- !duplicated(rid); last <- !duplicated(rid, fromLast = TRUE)
  sel_cls <- if (class == "all") rep(TRUE, length(rom)) else g$class_of_mut == class

  out <- do.call(rbind, lapply(g$samples, function(j) {
    sel <- sel_cls & g$sample_of_mut == j
    data.frame(region = seq_len(R), sample = j,
               observed = tabulate(rom[sel], nbins = R),
               expected = agg_regions(gof_lambda(g, class, j), rid),
               row.names = NULL, stringsAsFactors = FALSE)
  }))

  e <- pmax(out$expected, .Machine$double.eps)
  out$resid <- (out$observed - out$expected) / sqrt(e)
  out$excess <- out$observed - out$expected
  ## Half-counts, so a cell with zero observed or zero expected still has a
  ## finite fold change and stays on the volcano instead of being dropped.
  out$log2fc <- log2((out$observed + 0.5) / (out$expected + 0.5))
  if (compute_p) {
    out$p_val <- poisson_p_exact(out$observed, e)
    out$padj <- stats::p.adjust(out$p_val, method = p_adjust)
  }
  out$chrom <- g$chrom[first][out$region]
  out$start <- g$start[first][out$region]
  out$end <- g$end[last][out$region]
  out$class <- class
  out
}


#' Is a region's excess the cohort's, or one patient's?
#'
#' For every region, the share of the total POSITIVE excess that the single
#' largest contributor carries. Near 1 means the region is one patient's hotspot
#' and the pooled residual is not describing the cohort; near 1/nrow means the
#' excess is spread over the cohort and is a genuine failure of the intensity.
#'
#' Needed for both residual figures: a flagged region driven by one patient and a
#' flagged region the whole cohort shares are different findings and should not
#' be drawn the same way.
#'
#' @param rp Output of [gof_region_patient()].
#' @return A data frame with one row per region.
gof_region_driver <- function(rp) {
  exc <- pmax(rp$observed - rp$expected, 0)
  sp <- split(seq_len(nrow(rp)), rp$region)
  do.call(rbind, lapply(sp, function(i) {
    x <- exc[i]
    tot <- sum(x)
    k <- which.max(x)
    data.frame(region = rp$region[i][1],
               driver = rp$sample[i][k],
               driver_excess = x[k],
               driver_share = if (tot > 0) x[k] / tot else NA_real_,
               n_patients_pos = sum(x > 0),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}


#' Critical Kolmogorov-Smirnov distance, and what it means for power
#'
#' The two-sided asymptotic critical value, `D* = c(alpha)/sqrt(n)`. Quoted in
#' the QQ figure because it is the honest way to state what the equal-power test
#' can and cannot see: at `n` events the test rejects when the empirical
#' distribution of rescaled positions departs from the diagonal by more than
#' `D*` ANYWHERE. Power against a specific alternative is not computed and is not
#' needed - what the thinning buys is that `D*` is the SAME for every patient.
ks_critical <- function(n, alpha = 0.05) {
  c_alpha <- sqrt(-0.5 * log(alpha / 2))
  c_alpha / sqrt(n)
}


#' Rescaled-position QQ curves, on all of the mutations
#'
#' No thinning: every curve uses every mutation in its slice. Thinning belongs to
#' the burden comparison, not to a display of fit.
#'
#' TWO EQUIVALENT RENDERINGS are returned for each curve, and the equivalence is
#' exact but worth being precise about. Under the time-rescaling theorem the
#' rescaled positions are iid U(0, 1) given their number, and the gaps between
#' them are iid Exp(1); the gaps are a one-to-one transformation of the sorted
#' positions, so the two plots carry the SAME information and neither can show
#' something the other hides. What differs is what the eye picks up. The uniform
#' plot shows cumulative drift - a stretch of genome the model over- or
#' under-predicts pulls the whole curve off the diagonal and it stays off. The
#' exponential plot puts the deviations in the upper tail, where local clustering
#' (several mutations closer together than the intensity allows) shows as points
#' lifting off the line at the right-hand end.
#'
#' They are NOT interchangeable as tests: a KS statistic computed on the gaps is
#' a different statistic from one computed on the positions, and is the more
#' sensitive of the two against clustering. Only the plots are equivalent.
#'
#' @param g Output of [gof_setup()].
#' @param by "patient" (one curve each), "class" (one per macro class, patients
#'   pooled), "chromosome" (one per chromosome, patients pooled - the rescaling
#'   runs within the chromosome, so each curve is its own process), or "all"
#'   (the whole cohort as one superposed process).
#' @param n_points Points kept per curve. The extreme order statistics are always
#'   kept - the tail is where the departures are.
#' @param seed Set once for the randomised within-bin offsets.
#' @return A data frame with `group`, `n`, and both sets of coordinates.
gof_qq_curves <- function(g, by = c("patient", "class", "chromosome", "all"),
                          n_points = 400, seed = 1, min_events = 50,
                          verbose = TRUE) {
  by <- match.arg(by)
  set.seed(seed)
  groups <- switch(by, patient = g$samples, class = MUT_CLASSES,
                   chromosome = unique(g$chrom), all = "all")
  ## Computed once rather than per chromosome.
  lam_all <- if (by == "chromosome") gof_lambda(g, "all") else NULL
  chrom_of_mut <- if (by == "chromosome") g$chrom[g$bin_of_mut] else NULL

  out <- lapply(seq_along(groups), function(i) {
    grp <- groups[i]
    if (verbose && by == "patient" && i %% 20 == 0)
      message("  patient ", i, " / ", length(groups))
    if (by == "patient") {
      lam <- gof_lambda(g, "all", grp); ev <- g$bin_of_mut[g$sample_of_mut == grp]
    } else if (by == "class") {
      lam <- gof_lambda(g, grp); ev <- g$bin_of_mut[g$class_of_mut == grp]
    } else if (by == "chromosome") {
      ## Rescaled within the chromosome: Lambda accumulates over its bins only,
      ## so each curve is a self-contained Poisson process rather than a slice
      ## of one long concatenation.
      k <- which(g$chrom == grp)
      lam <- lam_all[k]
      ev <- match(g$bin_of_mut[chrom_of_mut == grp], k)
    } else {
      lam <- gof_lambda(g, "all"); ev <- g$bin_of_mut
    }
    ks <- gof_rescale_test(lam, ev, n0 = Inf, n_rep = 1, min_events = min_events)
    u <- ks$u_full
    n <- length(u)
    if (n < 2L) return(NULL)

    ## The rescaled interval between consecutive events, Delta_n = Lambda(x_n) -
    ## Lambda(x_{n-1}). Under the model these are iid Exp(1).
    delta <- diff(c(0, u)) * ks$lam_total

    ## KS PLOT, in the sense of Brown, Barbieri, Ventura, Kass & Frank (2002):
    ## the time-rescaling theorem says Delta_n ~ iid Exp(1), so
    ## z_n = 1 - exp(-Delta_n) ~ iid Uniform(0, 1); sorted and plotted against
    ## (n - 1/2)/N they should follow the 45-degree line, with Beta(n, N-n+1)
    ## bands from the order statistics.
    ##
    ## NOTE this is NOT the same transform as `unif_*` below. That one takes the
    ## rescaled POSITION Lambda(x_n)/Lambda(T), whose order statistics are
    ## uniform conditional on N. The KS plot works on the GAPS, which is the
    ## sharper of the two against local clustering, and is the named procedure.
    z <- sort(1 - exp(-delta))
    gaps <- sort(delta)
    ## Even-in-rank thinning straight-lines the upper tail: qexp explodes as the
    ## rank approaches n, so the final segment spans Exp(1) quantiles ~6 to ~14
    ## and draws a spurious elbow at qexp(1 - 1/n_points) in EVERY curve,
    ## whatever its size. The grid is refined geometrically at both ends.
    tail_ranks <- unique(round(exp(seq(0, log(max(n / 2, 2)),
                                       length.out = 150))))
    idx <- sort(unique(c(round(seq(1, n, length.out = min(n_points, n))),
                         tail_ranks, n + 1L - tail_ranks)))
    idx <- idx[idx >= 1 & idx <= n]
    data.frame(group = grp, n = n,
               unif_theoretical = (idx - 0.5) / n, unif_empirical = u[idx],
               exp_theoretical = stats::qexp((idx - 0.5) / n),
               exp_empirical = gaps[idx],
               ks_theoretical = (idx - 0.5) / n, ks_empirical = z[idx],
               idx = idx,
               ks_full = ks$ks_full, p_full = ks$p_full,
               ks_gap = ks$ks_gap, p_gap = ks$p_gap,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}


#' Draw the QQ curves from [gof_qq_curves()]
#'
#' @param df Output of [gof_qq_curves()].
#' @param scale "uniform" or "exponential"; see [gof_qq_curves()] on why both.
#' @param highlight Groups to colour and name. Defaults to the `n_highlight`
#'   worst by KS distance.
#' @param overlay A second [gof_qq_curves()] table drawn as one bold line on top
#'   - the cohort as a single superposed process, against which the per-patient
#'   curves are the spread.
#' @param overlay_colour Colour of that line.
#' @param band Draw the pointwise 95% Beta envelope of the order statistics. It
#'   is computed at the MEDIAN event count across the curves shown, so with
#'   curves of very different burden it is indicative rather than exact - at
#'   these sample sizes it is a hairline in any case, and a curve visibly off the
#'   diagonal is far outside it.
#' @param label_size,base_size Text sizes.
plot_gof_qq_curves <- function(df, scale = c("ks", "detrended", "uniform",
                                             "exponential"),
                               highlight = NULL, n_highlight = 5, band = TRUE,
                               only_highlight = FALSE, overlay = NULL,
                               overlay_colour = "darkblue",
                               overlay_label = "Total mutation",
                               label_size = 3, base_size = 11,
                               palette = NULL) {
  scale <- match.arg(scale)
  xv <- switch(scale, ks = , detrended = "ks_theoretical",
               uniform = "unif_theoretical", exponential = "exp_theoretical")
  yv <- switch(scale, ks = , detrended = "ks_empirical",
               uniform = "unif_empirical", exponential = "exp_empirical")
  df$.x <- df[[xv]]; df$.y <- df[[yv]]
  ## Detrended: the vertical distance from the 45-degree line rather than the
  ## quantile itself. On a [0, 1] axis these curves use ~17% of the panel height,
  ## so the departure is a barely visible bow; subtracting the diagonal spends
  ## the whole panel on it, and the KS distance becomes the largest excursion
  ## from zero, readable straight off the figure.
  if (scale == "detrended") df$.y <- df$.y - df$.x

  key <- unique(df[, c("group", "ks_full", "ks_gap", "n")])
  if (is.null(highlight))
    highlight <- utils::head(key$group[order(-key$ks_gap)], n_highlight)
  highlight <- highlight[!is.na(highlight)]
  hi <- df[df$group %in% highlight, , drop = FALSE]
  hi$group <- factor(hi$group, levels = highlight)
  rest <- if (only_highlight) df[0, , drop = FALSE]
          else df[!df$group %in% highlight, , drop = FALSE]

  p <- ggplot2::ggplot()
  if (band && nrow(key)) {
    n <- round(stats::median(key$n))
    k <- seq_len(min(n, 4000)); kk <- round(seq(1, n, length.out = length(k)))
    ## The band is the Beta(n, N-n+1) envelope of the order statistics. Its x
    ## coordinate follows the PANEL's x, which is the uniform quantile for both
    ## "ks" and "uniform" and only the exponential quantile for "exponential" -
    ## testing `scale == "uniform"` alone put exponential x on a uniform axis and
    ## drew a spurious curve from (0,0) to (1, 1 - 1/e).
    b <- data.frame(x = if (scale == "exponential") stats::qexp((kk - 0.5) / n)
                        else (kk - 0.5) / n,
                    lo = stats::qbeta(0.025, kk, n - kk + 1),
                    hi = stats::qbeta(0.975, kk, n - kk + 1))
    if (scale == "exponential") { b$lo <- -log(1 - b$lo); b$hi <- -log(1 - b$hi) }
    if (scale == "detrended") { b$lo <- b$lo - b$x; b$hi <- b$hi - b$x }
    p <- p + ggplot2::geom_ribbon(data = b,
                                  ggplot2::aes(x, ymin = lo, ymax = hi),
                                  fill = NA, outline.type = "full")
  }
  p <- p +
    (if (scale == "detrended")
       ggplot2::geom_hline(yintercept = 0, colour = "#CD2626", linewidth = 0.4)
     else ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#CD2626",
                               linewidth = 0.4)) +
    ggplot2::geom_line(data = rest, ggplot2::aes(.x, .y, group = group),
                       colour = "grey70", alpha = 0.5, linewidth = 0.5) +
    ggplot2::geom_line(data = hi, ggplot2::aes(.x, .y, colour = group),
                       linewidth = 0.55) +

    (if (is.null(overlay)) NULL else {
       overlay$.x <- overlay[[xv]]; overlay$.y <- overlay[[yv]]
       ## The overlay is built in its own branch, so the detrending has to be
       ## applied here too - otherwise it is drawn on the untransformed scale
       ## and runs diagonally across a panel whose curves are all near zero.
       if (scale == "detrended") overlay$.y <- overlay$.y - overlay$.x
       overlay$group <- factor(overlay_label, levels = c(highlight, overlay_label))
       ggplot2::geom_line(data = overlay, ggplot2::aes(.x, .y, colour = group),
                          linewidth = 0.55, linetype = "longdash")
     }) +
    ggplot2::labs(
      x = switch(scale, ks = , detrended = "Theoretical uniform quantiles",
                 uniform = "Uniform quantile",
                 exponential = "Exponential(1) quantile"),
      y = switch(scale,
                 detrended = "Observed - theoretical quantile",
                 ks = "Observed quantiles",
                 uniform = "Rescaled genomic position",
                 exponential = "Rescaled gap between mutations")) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
                   #legend.position = c(0.99, 0.01),
                   #legend.justification = c(1, 0),
                   #legend.background = ggplot2::element_rect(fill = "white",
                   #                                          colour = "grey70"),
                  # legend.key.size = ggplot2::unit(11, "pt"),
                   #legend.text = ggplot2::element_text(size = label_size * 2.6),
                   #legend.title = ggplot2::element_blank())
  ## Deliberately no blue among the patient colours: the pooled curve is light
  ## blue, and Set1's blue is close enough to it to be misread as a patient.
  if (is.null(palette))
    palette <- stats::setNames(
      rep(c("#9ECAE1", "#4FA3E3", "#2171B5", "#144C8C", "darkblue"),
          length.out = length(highlight)), highlight)
  if (!is.null(overlay)) palette <- c(palette, stats::setNames(overlay_colour,
                                                               overlay_label))
  p <- p + ggplot2::scale_colour_manual(
    values = palette, name = NULL,
    guide = if (length(highlight) || !is.null(overlay)) "legend" else "none")
  if (scale %in% c("ks", "uniform"))
    p <- p + ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1))

  p
}



#' Observed against model-implied distribution of per-bin counts
#'
#' The direct answer to "does the model capture the distribution of counts at
#' 2 kb": the ratio of observed to expected number of bins holding exactly m
#' mutations, computed in closed form rather than simulated.
#'
#' @param cd Output of [gof_count_distribution()], possibly rbind-ed over classes.
plot_gof_counts <- function(cd) {
  ## A class with no bin at all holding m mutations gives ratio 0, which a log
  ## axis cannot place. Those cells are the ones where the model expects a
  ## fraction of a bin and none is seen - nothing to read, so they are dropped
  ## rather than pinned to the axis.
  cd <- cd[is.finite(cd$ratio) & cd$ratio > 0, , drop = FALSE]
  cd$count <- factor(cd$count, levels = unique(cd$count))
  cd$class <- factor(cd$class, levels = intersect(c("all", MUT_CLASSES),
                                                  unique(cd$class)))
  ggplot2::ggplot(cd, ggplot2::aes(count, ratio, colour = class, group = class)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_colour_manual(values = MUT_COLS, name = NULL) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Mutations in a bin", y = "Observed / expected bins") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
}


#' Pearson residuals along the genome
#'
#' Chromosomes are drawn as alternating vertical bands, matching the
#' reconstruction track of Figure 3 so the two panels read as the same genome.
#' The percentage of megabase windows in each chromosome that depart from the
#' model is written along the top of the panel, above its band: a chromosome
#' with three huge residuals and one with a quarter of its windows off look
#' similar point by point and are different findings.
#'
#' What "depart" means: each window is tested on its own with the EXACT two-sided
#' Poisson tail probability of its count under the fitted intensity, and the
#' adjustment is Benjamini-Hochberg over every window in the genome. So the
#' per-chromosome percentages are shares of one genome-wide procedure, not 23
#' separate tests, and no chi-square is involved - `contrib`, the Pearson term
#' `(o - e)^2 / e`, is used only to rank regions and to build the effect sizes.
#'
#' Regions whose excess is carried by a SINGLE patient are drawn as open
#' circles. That distinction matters for reading the figure: a one-patient
#' hotspot is a statement about that patient, not about the intensity model, and
#' pooling has simply made it visible here.
#'
#' @param reg Output of [gof_regions()].
#' @param alpha Significance level for the flag.
#' @param n_label Regions to label, by absolute residual.
#' @param labels Optional character vector of labels, from [region_labels()],
#'   indexed by `reg$region`.
#' @param driver Optional output of [gof_region_driver()]. When supplied,
#'   flagged regions with `driver_share` above `driver_cut` are drawn open.
#' @param driver_cut Share of a region's positive excess above which it counts
#'   as one patient's.
#' @param pct_labels Write each chromosome's rejection rate along the top.
#' @param pct_size Text size for those labels.
#' @return A ggplot.
plot_gof_residual_track <- function(reg, alpha = 0.05, n_label = 10,
                                    labels = NULL, driver = NULL,
                                    driver_cut = 0.5, pct_labels = TRUE,
                                    pct_size = 3.4,
                                    band_cols = c("#4682B4", "antiquewhite"),
                                    point_cols = c("#000D8B", "#93BAF1")) {
  chr_lab <- unique(reg$chrom)
  num <- suppressWarnings(as.integer(sub("^chr", "", chr_lab)))
  chr_lab <- chr_lab[order(is.na(num), num, chr_lab)]
  reg$chr_i <- match(reg$chrom, chr_lab)
  reg <- reg[order(reg$chr_i, reg$start), ]
  reg$x <- seq_len(nrow(reg))
  reg$sig <- !is.na(reg$padj) & reg$padj < alpha
  reg$label <- if (is.null(labels)) "" else labels[reg$region]

  ## One patient's hotspot, or the cohort's?
  reg$one_patient <- FALSE
  if (!is.null(driver)) {
    d <- driver$driver_share[match(reg$region, driver$region)]
    reg$one_patient <- reg$sig & !is.na(d) & d > driver_cut & reg$resid > 0
  }

  ## Chromosome bands and their labels, from the region layout itself.
  bands <- do.call(rbind, lapply(seq_along(chr_lab), function(i) {
    k <- which(reg$chr_i == i)
    data.frame(chrom = chr_lab[i], xmin = min(reg$x[k]) - 0.5,
               xmax = max(reg$x[k]) + 0.5, mid = mean(range(reg$x[k])),
               n = length(k), n_sig = sum(reg$sig[k]),
               fill = factor(i %% 2), row.names = NULL,
               stringsAsFactors = FALSE)
  }))
  bands$pct <- 100 * bands$n_sig / bands$n
  xr <- c(0.5, nrow(reg) + 0.5)

  ## The x scale is applied ONCE, at the end: adding it here too makes ggplot
  ## replace it and warn on every panel.
  band_layer <- list(
    ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                                    ymax = Inf, fill = fill), alpha = 0.10),
    ggplot2::scale_fill_manual(values = band_cols, guide = "none"))

  ## ---- the track
  top <- reg[reg$sig, , drop = FALSE]
  top <- utils::head(top[order(-abs(top$resid)), , drop = FALSE], n_label)

  ## Headroom for the per-chromosome percentages, so they sit above every point
  ## rather than inside the cloud.
  rng <- range(reg$resid, na.rm = TRUE)
  ## Staggered on two rows: 23 chromosomes across a 9-inch panel leaves the
  ## narrow ones (19-22) closer together than their labels are wide, and a
  ## single row silently drops or overplots them.
  bands$row <- seq_len(nrow(bands)) %% 2
  y_pct <- rng[2] + c(0.09, 0.19)[bands$row + 1] * diff(rng)
  y_top <- rng[2] + 0.27 * diff(rng)

  p <- ggplot2::ggplot(reg, ggplot2::aes(x, resid)) +
    band_layer +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_point(data = reg[!reg$sig, ], colour = point_cols[2],
                        size = 0.4, alpha = 0.55) +
    ggplot2::geom_point(data = reg[reg$sig & !reg$one_patient, ],
                        colour = "#CD2626", size = 0.9, alpha = 0.95) +
    ggplot2::geom_point(data = reg[reg$one_patient, ], colour = "#CD2626",
                        size = 1.3, shape = 21, fill = NA, stroke = 0.5) +
    ggplot2::scale_x_continuous(limits = xr, expand = c(0.005, 0),
                                breaks = bands$mid,
                                labels = sub("^chr", "", chr_lab)) +
    ggplot2::labs(x = "Chromosome", y = "Pearson residual") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  if (pct_labels) {
    p <- p +
      ggplot2::geom_text(data = bands, inherit.aes = FALSE,
                         ggplot2::aes(x = mid, y = y_pct,
                                      label = sprintf("%.0f%%", pct)),
                         size = pct_size, colour = "black") +
      ggplot2::expand_limits(y = y_top)
  }

  if (nrow(top) && any(nzchar(top$label)) &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = top[nzchar(top$label), ], ggplot2::aes(label = label),
      size = 3.0, fontface = "italic", min.segment.length = 0,
      max.overlaps = 20, segment.colour = "grey45", segment.size = 0.25,
      box.padding = 0.25, ylim = c(NA, rng[2]))
  }
  p
}


#' Residual tracks, one facet per macro mutation class
#'
#' One ggplot faceted by class rather than six plots stacked, so the panels share
#' a genome axis and a residual scale by construction and the classes can be read
#' against each other. The strip carries the class and its rejection rate.
#'
#' The per-chromosome percentages are OFF here by default. At half panel width
#' twenty-three of them collide, and the strip already gives the rate that
#' matters for a class; the megabase-level percentages belong to the full-width
#' track in the main figure. `pct_labels = TRUE` turns them back on.
#'
#' @param reg_class Output of [gof_regions()] rbind-ed over `MUT_CLASSES`.
#' @param driver Optional named list of [gof_region_driver()] output, one entry
#'   per class; flagged regions one patient owns are then drawn open.
#' @param ncol Facet columns. The default 2 gives the 3 x 2 grid.
#' @param free_y Let each class set its own residual scale. Off by default: a
#'   shared scale is what makes "C>T is worse than T>G" visible rather than
#'   something the reader has to reconstruct from the axis numbers.
#' @inheritParams plot_gof_residual_track
plot_gof_residual_track_class <- function(reg_class, alpha = 0.05,
                                          labels = NULL, driver = NULL,
                                          driver_cut = 0.5, n_label = 3,
                                          classes = MUT_CLASSES, ncol = 2,
                                          pct_labels = FALSE, pct_size = 1.9,
                                          free_y = FALSE,
                                          band_cols = c("#4682B4", "antiquewhite"),
                                          point_cols = c("#000D8B", "#93BAF1")) {
  d <- reg_class[reg_class$class %in% classes, , drop = FALSE]

  ## One x ordering for every class: the regions are the same set, so the layout
  ## is computed once from one class and matched onto the rest.
  ref <- d[d$class == classes[1], , drop = FALSE]
  chr_lab <- unique(ref$chrom)
  num <- suppressWarnings(as.integer(sub("^chr", "", chr_lab)))
  chr_lab <- chr_lab[order(is.na(num), num, chr_lab)]
  ref$chr_i <- match(ref$chrom, chr_lab)
  ref <- ref[order(ref$chr_i, ref$start), ]
  xof <- stats::setNames(seq_len(nrow(ref)), ref$region)
  d$x <- unname(xof[as.character(d$region)])
  d$chr_i <- match(d$chrom, chr_lab)
  d$sig <- !is.na(d$padj) & d$padj < alpha

  d$one_patient <- FALSE
  if (!is.null(driver)) {
    for (cl in classes) {
      dv <- driver[[cl]]
      if (is.null(dv)) next
      k <- d$class == cl
      sh <- dv$driver_share[match(d$region[k], dv$region)]
      d$one_patient[k] <- d$sig[k] & d$resid[k] > 0 & !is.na(sh) & sh > driver_cut
    }
  }

  ## Strip text: class plus its rejection rate. Built as an ordered factor so the
  ## facets follow MUT_CLASSES and not alphabetical order.
  stat <- vapply(classes, function(cl) {
    k <- d$class == cl
    sprintf("%s   %d/%d Mb rejected (%.1f%%)", cl, sum(d$sig[k]), sum(k),
            100 * mean(d$sig[k]))
  }, character(1))
  d$facet <- factor(unname(stat[d$class]), levels = unname(stat))

  ref$x <- seq_len(nrow(ref))
  bands <- do.call(rbind, lapply(seq_along(chr_lab), function(i) {
    k <- which(ref$chr_i == i)
    data.frame(xmin = min(ref$x[k]) - 0.5, xmax = max(ref$x[k]) + 0.5,
               mid = mean(range(ref$x[k])), fill = factor(i %% 2),
               row.names = NULL)
  }))
  xr <- c(0.5, nrow(ref) + 0.5)

  p <- ggplot2::ggplot(d, ggplot2::aes(x, resid)) +
    ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                                    ymax = Inf, fill = fill), alpha = 0.10) +
    ggplot2::scale_fill_manual(values = band_cols, guide = "none") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_point(data = d[!d$sig, ], colour = point_cols[2],
                        size = 0.35, alpha = 0.55) +
    ggplot2::geom_point(data = d[d$sig & !d$one_patient, ], colour = "#CD2626",
                        size = 0.8, alpha = 0.95) +
    ggplot2::geom_point(data = d[d$one_patient, ], colour = "#CD2626",
                        size = 1.2, shape = 21, fill = NA, stroke = 0.5) +
    ggplot2::facet_wrap(~ facet, ncol = ncol,
                        scales = if (free_y) "free_y" else "fixed") +
    ggplot2::scale_x_continuous(limits = xr, expand = c(0.005, 0),
                                breaks = bands$mid,
                                labels = sub("^chr", "", chr_lab)) +
    ggplot2::labs(x = "Chromosome", y = "Pearson residual") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 8.5,
                                                      face = "bold"),
                   axis.text.x = ggplot2::element_text(size = 6))

  if (pct_labels) {
    pct <- do.call(rbind, lapply(classes, function(cl) {
      k <- d$class == cl
      dd <- d[k, ]
      do.call(rbind, lapply(seq_along(chr_lab), function(i) {
        j <- dd$chr_i == i
        data.frame(mid = bands$mid[i], pct = 100 * mean(dd$sig[j]),
                   row = i %% 2, facet = dd$facet[1], row.names = NULL)
      }))
    }))
    rng <- range(d$resid, na.rm = TRUE)
    pct$y <- rng[2] + c(0.09, 0.19)[pct$row + 1] * diff(rng)
    p <- p + ggplot2::geom_text(data = pct, inherit.aes = FALSE,
                                ggplot2::aes(x = mid, y = y,
                                             label = sprintf("%.0f%%", pct)),
                                size = pct_size, colour = "black") +
      ggplot2::expand_limits(y = rng[2] + 0.27 * diff(rng))
  }

  if (n_label > 0 && !is.null(labels) &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    top <- do.call(rbind, lapply(classes, function(cl) {
      dd <- d[d$class == cl & d$sig, , drop = FALSE]
      dd$label <- labels[dd$region]
      dd <- dd[nzchar(dd$label), , drop = FALSE]
      utils::head(dd[order(-abs(dd$resid)), , drop = FALSE], n_label)
    }))
    if (nrow(top))
      p <- p + ggrepel::geom_text_repel(
        data = top, ggplot2::aes(label = label), size = 2.4,
        fontface = "italic", min.segment.length = 0, max.overlaps = 20,
        segment.colour = "grey45", segment.size = 0.25, box.padding = 0.25)
  }

  ## Colour the strip text by class when ggh4x is available - the six colours are
  ## the cohort's mutation-class palette everywhere else in the paper. Wrapped in
  ## a guarded build: ggh4x's strip API has changed across versions, and a facet
  ## spec that errors at draw time would take the whole supplement with it.
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    alt <- try({
      q <- p + ggh4x::facet_wrap2(
        ~ facet, ncol = ncol, scales = if (free_y) "free_y" else "fixed",
        strip = ggh4x::strip_themed(
          text_x = ggh4x::elem_list_text(colour = unname(MUT_COLS[classes]),
                                         face = "bold", size = 8.5)))
      invisible(ggplot2::ggplot_build(q))
      q
    }, silent = TRUE)
    if (!inherits(alt, "try-error")) p <- alt
  }
  p
}


#' Equal-power uniform QQ plot, one line per patient
#'
#' Every patient is thinned to the SAME expected number of events before the
#' test, so the spread between the curves is a difference in fit and not a
#' difference in burden. `n0` must therefore be at or below the smallest burden
#' in the cohort: a patient with fewer events than `n0` cannot be thinned up to
#' it and is silently tested at lower power, which is the one way this figure can
#' lie. [gof_patients()] warns when that happens.
#'
#' The grey envelope is the pointwise 95% band that `n0` genuine uniforms would
#' produce, from the Beta order statistics. The dashed offsets are the
#' Kolmogorov-Smirnov critical distance at `n0` - the width of departure the test
#' can see, identical for every patient, which is the whole point of thinning.
#'
#' Only the worst few curves are coloured. Colouring every rejected patient
#' fills the panel and hides both the band and the bulk, and the bulk is the
#' finding: most patients track the diagonal.
#'
#' @param qq The `qq` element of [gof_patients()]`(qq = TRUE)`.
#' @param pat The `patients` element, supplying the ranking and the test.
#' @param alpha Significance level.
#' @param n_highlight Curves to colour and name, worst first by KS distance.
plot_gof_qq <- function(qq, pat, alpha = 0.05, n_highlight = 5,
                        rank_by = "focal") {
  n <- round(stats::median(pat$n_thin, na.rm = TRUE))
  worst <- utils::head(pat$sample[order(-pat[[rank_by]])], n_highlight)
  df <- qq[qq$sample %in% worst, , drop = FALSE]
  df$sample <- factor(df$sample, levels = worst)
  rest <- qq[!qq$sample %in% worst, , drop = FALSE]

  k <- seq_len(n)
  band <- data.frame(theoretical = (k - 0.5) / n,
                     lo = stats::qbeta(0.025, k, n - k + 1),
                     hi = stats::qbeta(0.975, k, n - k + 1))
  Dstar <- ks_critical(n, alpha)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band,
                         ggplot2::aes(theoretical, ymin = lo, ymax = hi),
                         fill = "grey85") +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey45",
                         linewidth = 0.3) +
    ggplot2::geom_abline(slope = 1, intercept = c(-Dstar, Dstar), linetype = 2,
                         colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_line(data = rest,
                       ggplot2::aes(theoretical, empirical, group = sample),
                       colour = "grey40", alpha = 0.22, linewidth = 0.25) +
    ggplot2::geom_line(data = df,
                       ggplot2::aes(theoretical, empirical, colour = sample),
                       linewidth = 0.6) +
    ggplot2::scale_colour_brewer(palette = "Set1", name = NULL) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = "Uniform quantile", y = "Rescaled genomic position") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = c(0.99, 0.01),
                   legend.justification = c(1, 0),
                   legend.background = ggplot2::element_rect(fill = "white",
                                                             colour = "grey70"),
                   legend.key.size = ggplot2::unit(11, "pt"),
                   legend.text = ggplot2::element_text(size = 8.5))
}


#' Volcano of (patient, region) cells
#'
#' Which patient, in which megabase, the model gets most wrong - the resolution
#' the pooled residual track cannot reach, because a region every patient is
#' slightly off in and a region one patient is wildly off in look the same once
#' the counts are summed.
#'
#' The x axis is observed MINUS expected on a signed log modulus scale
#' (John & Draper 1980), `sign(x) * log10(1 + |x|)`, ticked in real mutation
#' counts. It ranks by how many mutations are actually unexplained rather than by
#' fold change: 95% of cells sit within 10 mutations of zero while the tail
#' reaches +357, a linear axis collapses that into a spike, and a plain log axis
#' cannot carry zero or the sign.
#'
#' There is NO expected-count filter. Filtering on expected count is not neutral
#' here - a large excess in a low-expectation cell is what a hotspot is, so a
#' threshold deletes the largest departures first. An earlier `min_expected = 8`,
#' set for label spacing, removed DO1020's 363-against-5.5 kataegis from the
#' figure meant to show it.
#'
#' Points are labelled `PATIENT` over `chrN:M`, M being the megabase within that
#' chromosome. No gene symbols: naming the longest gene in a megabase implies a
#' precision the window does not have.
#'
#' @param rp Output of [gof_region_patient()].
#' @param alpha Significance level, applied to `rp$padj`.
#' @param n_label Cells to name, split between the two arms.
#' @param y_cap `-log10 p` above this is truncated to it and drawn as triangles,
#'   with a rule at the cut. A few cells reach 300 and flatten everything else
#'   into the bottom eighth of the panel.
#' @param label_floor,point_size Layout.
#' @param headroom Panel top, as a multiple of `y_cap`. The space above the rule
#'   is where the labels of the capped points go.
plot_gof_volcano <- function(rp, alpha = 0.05, n_label = 12, y_cap = 100,
                             point_size = 0.6, label_floor = 0.30,
                             headroom = 1.34) {
  d <- rp
  d$sig <- !is.na(d$padj) & d$padj < alpha
  ## The SIGN is the x axis - colouring by it as well says the same thing twice.
  ## Colour carries the one thing position does not: whether the cell is rejected.
  d$dir <- factor(ifelse(!d$sig, "not rejected", "rejected"),
                  levels = c("not rejected", "rejected"))
  d$arm <- ifelse(d$log2fc > 0, "excess", "deficit")   # label placement only
  d$logp <- -log10(pmax(d$p_val, .Machine$double.xmin))
  d$tag <- paste0(d$sample, "\n", d$chrom, ":", round(d$start / 1e6))

  cap_used <- !is.null(y_cap) && any(d$logp > y_cap)
  ## Headroom above the cap line. The capped points all sit at exactly y_cap, so
  ## their labels have nowhere to go inside the panel and collide with whatever
  ## else is near the top; opening space above the rule gives repel somewhere to
  ## put them. `headroom` is the factor of y_cap the panel extends to.
  y_top <- if (cap_used) y_cap * headroom else NA_real_
  d$clipped <- if (cap_used) d$logp > y_cap else FALSE
  d$logp_plot <- if (cap_used) pmin(d$logp, y_cap) else d$logp

  slog <- function(x) sign(x) * log10(1 + abs(x))
  d$.x <- slog(d$excess); d$.y <- d$logp_plot

  lab <- do.call(rbind, lapply(c("excess", "deficit"), function(k) {
    sel <- d$sig & d$arm == k
    dk <- d[sel, , drop = FALSE]
    if (!nrow(dk)) return(NULL)
    utils::head(dk[order(-abs(d$excess[sel])), ], ceiling(n_label / 2))
  }))

  p <- ggplot2::ggplot(d, ggplot2::aes(.x, .y, colour = dir)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3)
  if (cap_used)
    p <- p + ggplot2::geom_hline(yintercept = y_cap, linetype = 2,
                                 colour = "grey40", linewidth = 0.35) +
      ggplot2::expand_limits(y = y_top)

  br <- c(-300, -100, -30, -10, -3, 0, 3, 10, 30, 100, 300)
  br <- br[abs(br) <= max(abs(d$excess)) * 1.5]

  p <- p +
    ggplot2::geom_point(data = d[!d$clipped, ], size = point_size, alpha = 0.45,
                        shape = 16) +
    ggplot2::geom_point(data = d[d$clipped, ], size = point_size * 2.6,
                        shape = 17) +
    ggplot2::scale_colour_manual(
      values = c("not rejected" = "#C9D6E8", "rejected" = "#CD2626"),
      name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = slog(br),
      labels = ifelse(br > 0, paste0("+", br), as.character(br))) +
    ## Ticks stop at the cap: the headroom above it exists for labels, and a
    ## tick at 120 would imply data that was truncated away.
    (if (cap_used)
       ggplot2::scale_y_continuous(breaks = pretty(c(0, y_cap)))
     else NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 2.2, alpha = 1))) +
    ggplot2::labs(x = "Observed - expected mutation each patient (1Mb)",
                  y = expression(-log[10](p))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   aspect.ratio = 1,
                   legend.position = "none")

  if (!is.null(lab) && nrow(lab) && requireNamespace("ggrepel", quietly = TRUE)) {
    yhi <- max(d$logp_plot, na.rm = TRUE)
    xr <- range(d$.x, na.rm = TRUE)
    ## force/box.padding matter more than max.overlaps: a high overlap budget
    ## makes ggrepel DRAW colliding labels rather than drop them.
    common <- list(size = 2.9, colour = "black", lineheight = 0.85,
                   min.segment.length = 0, max.overlaps = 40,
                   segment.size = 0.2, segment.colour = "grey55",
                   box.padding = 0.7, point.padding = 0.3,
                   force = 12, force_pull = 0.6, show.legend = FALSE)
    left <- lab[lab$arm == "deficit", , drop = FALSE]
    right <- lab[lab$arm == "excess", , drop = FALSE]
    ## The deficit arm stops below the legend, yielding the one clean corner.
    if (nrow(left))
      p <- p + do.call(ggrepel::geom_text_repel, c(list(
        data = left, mapping = ggplot2::aes(label = tag),
        xlim = c(xr[1], -0.05), ylim = c(label_floor * yhi, 0.82 * yhi),
        direction = "both"), common))
    if (nrow(right))
      p <- p + do.call(ggrepel::geom_text_repel, c(list(
        data = right, mapping = ggplot2::aes(label = tag),
        xlim = c(0.05, xr[2]),
        ylim = c(label_floor * yhi, if (cap_used) y_top else NA),
        direction = "both"), common))
  }
  ## The patients named here, so a companion panel can highlight exactly these
  ## and not a different set chosen by a second, silently diverging rule.
  attr(p, "labelled") <- if (is.null(lab)) character(0) else unique(lab$sample)
  p
}
