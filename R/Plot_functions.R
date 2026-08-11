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
