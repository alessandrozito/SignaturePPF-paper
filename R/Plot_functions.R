## Figure helpers shared across analyses.

#' Colours for the signatures reported in the breast analyses, grouped by aetiology
SIG_COLS <- c(
  SBS1   = "#E31A1C",   # clock
  SBS5   = "#FDBF6F",   # clock
  SBS2   = "#6A3D9A",   # APOBEC
  SBS13  = "#111111",   # APOBEC
  SBS3   = "#1F78B4",   # HRD
  SBS6   = "#33A02C",   # MMRd
  SBS20  = "#7FBC41",   # MMRd
  SBS26  = "#B2DF8A",   # MMRd
  SBS44  = "#00694E",   # MMRd
  SBS8   = "#00CFE0",   # ROS
  SBS18  = "#0000FF",   # ROS
  SBS17a = "#FB9A99",   # unknown
  SBS17b = "#E7298A",   # unknown
  SBS30  = "#FFC0CB",   # BER
  SBS40a = "#FF00FF"    # unknown / flat
)


#' Colours for a set of signatures, falling back for anything not in SIG_COLS
#'
#' SIG_COLS only names the COSMIC signatures the breast analyses report. A de
#' novo fit labels its signatures SigN01, SigN02, ... and would otherwise fall
#' through the manual scale with a warning and no colour at all.
sig_palette <- function(sigs) {
  out <- SIG_COLS[sigs]
  gap <- is.na(out)
  if (any(gap)) out[gap] <- grDevices::hcl.colors(sum(gap), "Dark 3")
  stats::setNames(unname(out), sigs)
}


#' One panel per signature: relevance weight against mutations attributed
#'
#' A compact companion to a beta heatmap. Point size is \eqn{\mu_k}, fill is the
#' number of mutations the signature wins, and a cross marks signatures the
#' compressive prior has switched off, so a coefficient row can be read together
#' with whether its signature is actually carrying anything.
#'
#' @param df_Assign A data frame from [df_assign()].
#' @param levs Signature order, normally `colnames(fit$Betas)` so the panels line
#'   up row-for-row with a beta heatmap plotted beside it.
#' @param mu_cutoff Below this `mu`, a signature is drawn as compressed.
plot_vector_facets_x <- function(df_Assign, levs, mu_cutoff = 0.05) {
  df_Assign$best_sig <- factor(df_Assign$best_sig, levels = levs)
  df_Assign$compressed <- ifelse(df_Assign$m == 0 | df_Assign$mu < mu_cutoff,
                                 "compressed", "not compressed")

  ggplot2::ggplot(df_Assign, ggplot2::aes(x = 1, y = 1, size = .data$mu,
                                   fill = .data$m, shape = .data$compressed)) +
    ggplot2::geom_point(colour = "black", stroke = 0.7) +
    ggplot2::facet_wrap(~ best_sig, ncol = 1, strip.position = "left") +
    ggplot2::scale_size(name = expression(mu[k]), range = c(3, 12)) +
    ggplot2::scale_fill_gradientn(
      name = "N. mutations",
      colours = c("#F6C866", "#F2AB67", "#EF8F6B", "#ED7470", "#BF6E97",
                  "#926AC2", "#6667EE", "#4959C7", "#2D4A9F", "#173C78")) +
    ggplot2::scale_shape_manual(
      name = "Compressed",
      values = c("compressed" = 4, "not compressed" = 21)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      strip.text.y.left = ggplot2::element_text(angle = 0),
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "left",
      panel.spacing = grid::unit(0.03, "lines"))
}


#' Coefficient agreement between two cohorts
#'
#' One point per (signature, covariate), the coefficient in cohort A against
#' cohort B. Signatures switched off in BOTH cohorts are dropped - their
#' coefficients are prior draws and comparing them measures nothing - and
#' signatures present in only one are drawn hollow, since they are not
#' comparable either.
#'
#' @param fit_x,fit_y The two fits.
#' @param label_x,label_y Axis labels.
#' @param mu_tol A signature counts as present when `mu` exceeds this.
#' @param inset Add a zoomed inset over the region where most coefficients lie.
plot_beta_replication <- function(fit_x, fit_y, label_x = "Cohort A",
                                  label_y = "Cohort B", mu_tol = 0.01,
                                  inset = TRUE, zoom = 0.1) {
  sigs <- intersect(colnames(fit_x$Betas), colnames(fit_y$Betas))
  present_x <- fit_x$Mu[sigs] > mu_tol
  present_y <- fit_y$Mu[sigs] > mu_tol
  keep <- sigs[present_x | present_y]
  if (!length(keep)) stop("no signature is present in either cohort")

  d <- do.call(rbind, lapply(keep, function(s) {
    data.frame(sig = s,
               feature = rownames(fit_x$Betas),
               x = fit_x$Betas[, s],
               y = fit_y$Betas[, s],
               stringsAsFactors = FALSE)
  }))
  shared <- keep[present_x[keep] & present_y[keep]]
  d$status <- ifelse(!(d$sig %in% shared), "present in one cohort",
                     ifelse(sign(d$x) == sign(d$y), "sign agrees", "sign flips"))
  d$sig <- factor(d$sig, levels = keep)

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linetype = 3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55", linetype = 3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$sig, shape = .data$status),
                        size = 1, stroke = 1) +
    ggplot2::scale_shape_manual(values = c("sign agrees" = 19,
                                           "sign flips" = 4,
                                           "present in one cohort" = 1)) +
    ggplot2::scale_colour_manual(values = sig_palette(levels(d$sig)),
                                 breaks = levels(d$sig)) +
    ggplot2::labs(x = bquote(beta ~ "(" * .(label_x) * ")"),
                  y = bquote(beta ~ "(" * .(label_y) * ")"),
                  colour = "Signature", shape = NULL) +
    ggplot2::theme_bw() +
    ggplot2::annotate("rect", xmin = -zoom, xmax = zoom, ymin = -zoom, ymax = zoom,
                      fill = NA, colour = "grey40", linetype = 2)

  if (!inset) return(p)

  p_zoom <- p +
    ggplot2::coord_fixed(xlim = c(-zoom, zoom), ylim = c(-zoom, zoom)) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(legend.position = "none",
                   axis.title = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white",
                                                           colour = "grey40"),
                   plot.margin = ggplot2::margin(2, 2, 2, 2))
  p + patchwork::inset_element(p_zoom, left = 0.55, bottom = 0.02,
                               right = 0.98, top = 0.45, align_to = "panel")
}


#' Coefficient of each covariate along a sequence of nested models
#'
#' One panel per covariate, showing what happens to its coefficient as further
#' covariates are added to the model. A coefficient that moves little across the
#' sequence is one the other covariates do not explain away.
#'
#' Only signatures the compressive prior keeps in EVERY model are drawn. For a
#' switched-off signature \eqn{\beta} is a draw from its prior, so a trajectory
#' that includes one would show the prior wandering rather than an estimate
#' changing.
#'
#' @param fits The fitted models, in the order the covariates were added. Each
#'   must be a fit whose `Betas` rows are the covariates it was given.
#' @param covariate_order The covariate added at each step, in order.
#' @param mu_cutoff A signature counts as present when `mu` exceeds this.
plot_beta_path <- function(fits, covariate_order, mu_cutoff = 0.05) {
  d <- do.call(rbind, lapply(seq_along(fits), function(m) {
    B <- fits[[m]]$Betas
    data.frame(model = m,
               covariate = rep(rownames(B), times = ncol(B)),
               signature = rep(colnames(B), each = nrow(B)),
               beta = as.numeric(B),
               mu = rep(as.numeric(fits[[m]]$Mu[colnames(B)]), each = nrow(B)),
               stringsAsFactors = FALSE)
  }))

  alive <- vapply(split(d$mu, d$signature), function(z) all(z > mu_cutoff),
                  logical(1))
  d <- d[d$signature %in% names(which(alive)), , drop = FALSE]
  if (!nrow(d)) stop("no signature is present in every model")
  d$covariate <- factor(d$covariate, levels = covariate_order)
  d$signature <- factor(d$signature, levels = sort(unique(d$signature)))

  ggplot2::ggplot(d, ggplot2::aes(.data$model, .data$beta,
                                  colour = .data$signature,
                                  group = .data$signature)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::facet_grid(~ covariate, scales = "free", space = "free_x") +
    ggplot2::scale_colour_manual(values = sig_palette(levels(d$signature))) +
    ggplot2::scale_y_continuous(n.breaks = 8) +
    ggplot2::scale_x_continuous(breaks = seq_along(covariate_order),
                                labels = paste0("+", covariate_order),
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::labs(x = "Model (covariates added sequentially)",
                  y = expression(hat(beta)), colour = "Signature") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(angle = 90),
                   strip.clip = "off",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                       size = 7))
}


#' How mutation attribution flows as covariates are added
#'
#' An alluvial diagram over a sequence of models: each stratum is the set of
#' mutations a model attributes to one signature, and a ribbon between two
#' columns is a set of mutations that moved. The percentage above each column is
#' the share of mutations that changed signature when that covariate entered.
#'
#' @param A A character matrix, mutations in rows and models in columns, holding
#'   the signature each model attributes each mutation to.
#' @param labels Column labels, length `ncol(A)`.
#' @param levs Signature order for the strata and the legend.
plot_assignment_alluvial <- function(A, labels, levs = NULL) {
  if (is.null(levs)) levs <- sort(unique(as.vector(A)))
  n_models <- ncol(A)

  # Mutations sharing a whole path through the sequence are interchangeable, so
  # they are collapsed to one ribbon carrying a count. Without this the diagram
  # would be drawn from several hundred thousand rows.
  key <- do.call(paste, c(as.data.frame(A, stringsAsFactors = FALSE),
                          list(sep = "\r")))
  tab <- table(key)
  paths <- do.call(rbind, strsplit(names(tab), "\r", fixed = TRUE))

  lodes <- data.frame(
    path_id = rep(seq_len(nrow(paths)), times = n_models),
    x = rep(seq_len(n_models), each = nrow(paths)),
    stratum = factor(as.vector(paths), levels = levs),
    n = rep(as.integer(tab), times = n_models))

  moved <- vapply(seq_len(n_models - 1L),
                  function(j) mean(A[, j] != A[, j + 1L]), numeric(1))
  lab_df <- data.frame(x = seq_len(n_models - 1L) + 1L, y = nrow(A),
                       lab = sprintf("%.1f%%", 100 * moved))

  ggplot2::ggplot(lodes,
                  ggplot2::aes(x = .data$x, stratum = .data$stratum,
                               alluvium = .data$path_id, y = .data$n,
                               fill = .data$stratum)) +
    ggalluvial::geom_flow(alpha = 0.6, width = 0.3) +
    ggalluvial::geom_stratum(width = 0.3, colour = "grey30", linewidth = 0.2) +
    ggplot2::geom_text(data = lab_df, inherit.aes = FALSE,
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$lab),
                       vjust = -0.4, size = 3) +
    ggplot2::scale_fill_manual(values = sig_palette(levs), name = "Signature",
                               drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_models), labels = labels) +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(x = "Model (covariates added sequentially)", y = "Mutations") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Per-patient predictive error along a sequence of nested models
#'
#' One box per model and sample split, so the in-sample and out-of-sample curves
#' can be read against each other: covariates that only fit noise improve the
#' first while leaving the second flat or worse.
#'
#' @param rmse A data frame with `model`, `patient`, `in_sample`, `out_sample`.
#' @param labels Model labels, one per level of `model`.
plot_rmse_path <- function(rmse, labels) {
  d <- data.frame(
    model = factor(rep(rmse$model, 2), levels = sort(unique(rmse$model))),
    patient = rep(rmse$patient, 2),
    rmse = c(rmse$in_sample, rmse$out_sample),
    set = rep(c("In-sample", "Out-of-sample"), each = nrow(rmse)))

  ggplot2::ggplot(d, ggplot2::aes(.data$model, .data$rmse, colour = .data$set)) +
    ggplot2::geom_boxplot(outlier.shape = NA,
                          position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(
      jitter.width = 0.15, dodge.width = 0.8), size = 0.7, alpha = 0.4) +
    ggplot2::scale_x_discrete(labels = labels) +
    ggplot2::scale_colour_manual(values = c("In-sample" = "#4682B4",
                                            "Out-of-sample" = "#CD2626")) +
    ggplot2::labs(x = "Model", y = "Per-patient RMSE (1 Mb regions)",
                  colour = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")
}


#' Mutation burden per megabase along the genome, for one or more cohorts
#'
#' @param datasets A named list of cohort objects.
#' @param window Window width in bases.
plot_burden_along_genome <- function(datasets, window = 1e6) {
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  bins <- GenomicRanges::tileGenome(
    GenomeInfoDb::seqlengths(genome)[paste0("chr", c(1:22, "X"))],
    tilewidth = window, cut.last.tile.in.chrom = TRUE)

  df <- do.call(rbind, lapply(names(datasets), function(nm) {
    gr <- datasets[[nm]]$gr_Mutations
    region <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(gr, bins))
    tab <- table(factor(region, levels = seq_along(bins)))
    data.frame(region = seq_along(bins), n = as.integer(tab), data = nm)
  }))

  bands <- as.data.frame(bins)
  bands$region <- seq_along(bins)
  bands <- do.call(rbind, lapply(split(bands, bands$seqnames), function(z) {
    data.frame(seqnames = z$seqnames[1], xmin = min(z$region),
               xmax = max(z$region) + 1)
  }))
  bands <- bands[order(bands$xmin), ]
  bands$fill_color <- factor(seq_len(nrow(bands)) %% 2)

  ggplot2::ggplot() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::geom_rect(data = bands,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = 0, ymax = Inf, fill = .data$fill_color),
                       alpha = 0.2) +
    ggplot2::scale_fill_manual(values = c("#4682B4", "antiquewhite2"),
                               guide = "none") +
    ggplot2::geom_line(data = df,
                       ggplot2::aes(.data$region, .data$n, colour = .data$data),
                       linewidth = 0.3) +
    ggplot2::labs(x = "Genomic position (Mb bins)", y = "Mutations",
                  colour = NULL) +
    ggplot2::theme(legend.position = "bottom")
}
