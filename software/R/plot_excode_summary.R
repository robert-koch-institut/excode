#' Plot EXCODE summary diagnostics across multiple panels
#'
#' @description
#' Visualizes an EXCODE summary as a multi-panel figure showing
#' (i) observed counts with model fits (`mu*`),
#' (ii) excess probability \eqn{1 - \mathrm{posterior0}},
#' (iii) z-scores,
#' (iv) \eqn{-\log_{10}(p\text{-value})}, and
#' (v) threshold “tiles” indicating which metrics exceed user-defined cutoffs.
#'
#' @param excode_summary A data frame (or tibble) with at least the columns:
#'   \itemize{
#'     \item \code{date}: Date or POSIXt.
#'     \item \code{observed}: Numeric observed counts.
#'     \item \code{posterior0}: Numeric posterior probability of \emph{no} excess (can contain \code{NA}, treated as 1).
#'     \item \code{zscore}: Numeric z-score (can contain \code{NA}, treated as 0).
#'     \item \code{pval}: Numeric p-value (can contain \code{NA}, treated as 1).
#'     \item One or more columns named \code{mu0}, \code{mu1}, … with model fits.
#'   }
#' @param states Optional integer giving the number of latent states (i.e., how
#'   many \code{mu*} series to consider). If missing, it is inferred from the
#'   available columns matching \code{"^mu\\d+$"}. If provided, only
#'   \code{mu0} through \code{mu{states-1}} that exist in \code{excode_summary}
#'   are used.
#' @param posterior_cutoff Numeric in \eqn{[0,1]}. Horizontal reference line and
#'   threshold for the excess probability panel (default \code{0.5}).
#' @param zscore_cutoff Numeric. Horizontal reference line and threshold for the
#'   z-score panel (default \code{2}).
#' @param p_cutoff Numeric in \eqn{(0,1]}. Transformed to \eqn{-\log_{10}}
#'   for the p-value panel and used as a threshold in the tiles (default \code{0.025}).
#' @param type Character, either \code{"line"} or \code{"bar"} (partial matching
#'   disabled via \code{match.arg}). Controls whether panels are drawn with
#'   ribbons/lines (\code{"line"}) or columns (\code{"bar"}). In the first panel,
#'   \code{"bar"} draws \code{observed} as bars and \code{mu*} as steps; \code{"line"}
#'   draws all series as lines.
#'
#' @details
#' Missing values are handled defensively before plotting:
#' \code{posterior0 = 1} (no excess), \code{zscore = 0}, \code{pval = 1}.
#' A consistent black/green/orange-red palette is used for the first panel:
#' observed in black; \code{mu0} in green; remaining \code{mu*} from gold to tomato.
#' Facet strips are placed on the right and a common theme is applied so panels
#' align cleanly when combined with \pkg{patchwork}.
#'
#' The final plot stacks five panels (from top to bottom):
#' \enumerate{
#'   \item Observed counts and model fits.
#'   \item Excess probability \eqn{1 - \mathrm{posterior0}} with a dotted cutoff.
#'   \item z-score with a dotted cutoff (ribbon handles positive/negative).
#'   \item \eqn{-\log_{10}(p)} with a dotted cutoff at \eqn{-\log_{10}(\mathrm{p\_cutoff})}.
#'   \item Binary “tiles” summarizing whether each metric exceeds its cutoff.
#' }
#'
#' @return A \code{patchwork} object (composed of \pkg{ggplot2} plots) that can be
#' further modified with \code{+} (e.g., themes, scales) or saved with \code{ggsave()}.
#'
#' @section Required packages:
#' Relies on \pkg{ggplot2}, \pkg{dplyr}, \pkg{tidyr}, \pkg{patchwork},
#' \pkg{grid}, and \pkg{grDevices}.
#'
#' @seealso
#' \code{\link[patchwork]{plot_layout}}, \code{\link[ggplot2]{ggplot}},
#' \code{\link[dplyr]{mutate}}, \code{\link[tidyr]{pivot_longer}}
#' 
#' @examples
#'
#' data(mort_df_germany)
#' res_har_nb <- run_excode(surv_ts = mort_df_germany,
#' timepoints = 325,
#' distribution = "NegBinom",
#' states = 3,
#' periodic_model = "Harmonic",
#' time_trend = "Spline2",
#' return_full_model = TRUE) 
#' 
#' sum_har_nb <- summary(res_har_nb)
#' 
#' plot_excode_summary(sum_har_nb, type="line")
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_col geom_step geom_ribbon geom_hline
#'   scale_color_manual facet_wrap xlab ylab theme_bw theme element_blank guides guide_legend
#' @importFrom dplyr select rename mutate bind_rows case_when all_of
#' @importFrom tidyr pivot_longer
#' @importFrom patchwork plot_layout
#' @importFrom grid unit
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom rlang .data .env
#' @export
plot_excode_summary <- function(
  excode_summary,
  states,
  posterior_cutoff = 0.5,
  zscore_cutoff = 2,
  p_cutoff = 0.025,
  type = c("line", "bar")
) {
  type <- match.arg(type)

  # clean NAs safely (stay in data mask; no .data needed)
  excode_summary <- excode_summary %>%
    dplyr::mutate(
      posterior0 = dplyr::coalesce(excode_summary$posterior0, 1),
      zscore     = dplyr::coalesce(excode_summary$zscore, 0),
      pval       = dplyr::coalesce(excode_summary$pval, 1)
    )

  # infer mu columns / states
  mu_names <- grep("^mu\\d+$", names(excode_summary), value = TRUE)
  if (missing(states)) {
    states <- length(mu_names)
  } else {
    mu_names <- paste0("mu", 0:(states - 1))
    mu_names <- intersect(mu_names, names(excode_summary))
  }

  # Common theme to keep appearance identical and ensure right-side strips show
  panel_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      strip.placement = "outside",
      plot.margin = grid::unit(c(0, 0.5, 0, 0), "cm")
    )

  # -------------------- Panel 1: Observed + mu fits --------------------
  df_obs <- excode_summary %>%
    dplyr::select(date, observed) %>%
    dplyr::rename(value = observed) %>%
    dplyr::mutate(series = "Observed")

  df_mus <- excode_summary %>%
    dplyr::select(date, dplyr::all_of(mu_names)) %>%
    tidyr::pivot_longer(-date, names_to = "series", values_to = "value")

  df_p1 <- dplyr::bind_rows(df_obs, df_mus) %>%
    dplyr::mutate(group = "Number of cases\nand model fit")

  # Color palette: black observed; greens for mu0, gold→tomato for remaining mu
  mu_cols <- c(
    grDevices::colorRampPalette(c("lightgreen", "forestgreen"))(5)[2],
    if (length(mu_names) > 1) grDevices::colorRampPalette(c("gold", "tomato"))(length(mu_names) - 1) else character(0)
  )
  names(mu_cols) <- mu_names
  legend_values <- c("Observed" = "black", mu_cols)

  p1 <- ggplot2::ggplot(df_p1, ggplot2::aes(x = .data$date, y = .data$value, color = .data$series))
  if (type == "bar") {
    # Observed as bars; mu* as steps
    p1 <- p1 +
      ggplot2::geom_col(
        data = subset(df_p1, df_p1$series == "Observed"),
        fill = "lightgrey", color = "grey", alpha = 0.9
      ) +
      ggplot2::geom_step(
        data = subset(df_p1, df_p1$series != "Observed"),
        linewidth = 0.9
      )
  } else {
    # All as lines (mu* will be on top)
    p1 <- p1 + ggplot2::geom_line(linewidth = 0.9)
  }
  p1 <- p1 +
    ggplot2::scale_color_manual(values = legend_values, breaks = names(legend_values)) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    panel_theme +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.key.width = grid::unit(1, "cm")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, byrow = TRUE))

  # -------------------- Panel 2: Posterior (1 - posterior0) --------------------
  df2 <- excode_summary %>%
    dplyr::mutate(
      excess_prob = 1 - excode_summary$posterior0, # bare names are fine inside mutate
      group       = "Excess\nProbability"
    )

  if (type == "line") {
    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = .data$date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$excess_prob),
        fill = "lightblue", color = "lightblue", alpha = 0.75
      )
  } else {
    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = .data$date, y = .data$excess_prob)) +
      ggplot2::geom_col(fill = "lightblue", color = "lightblue", alpha = 0.75)
  }
  p2 <- p2 +
    ggplot2::geom_hline(yintercept = posterior_cutoff, linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::scale_y_continuous(breaks = c(0, posterior_cutoff, 1)) +
    panel_theme

  # -------------------- Panel 3: zscores --------------------
  # For ribbons, handle negatives so ribbons render both above and below zero
  df3 <- excode_summary %>%
    dplyr::mutate(
      ymin_zscore = pmin(0, excode_summary$zscore),
      ymax_zscore = pmax(0, excode_summary$zscore),
      group       = "zscore\n"
    )

  if (type == "line") {
    p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = .data$date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ymin_zscore, ymax = .data$ymax_zscore),
        fill = "cornflowerblue", color = "cornflowerblue", alpha = 0.75
      )
  } else {
    p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = .data$date, y = .data$zscore)) +
      ggplot2::geom_col(fill = "cornflowerblue", color = "cornflowerblue", alpha = 0.75)
  }
  p3 <- p3 +
    ggplot2::geom_hline(yintercept = zscore_cutoff, linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    panel_theme

  # -------------------- Panel 4: -log10(pval) --------------------
  df4 <- excode_summary %>%
    dplyr::mutate(
      logp = -log10(excode_summary$pval),
      group = "-log10(pval)\n"
    )

  if (type == "line") {
    p4 <- ggplot2::ggplot(df4, ggplot2::aes(x = .data$date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$logp),
        fill = "darkblue", color = "darkblue", alpha = 0.75
      )
  } else {
    p4 <- ggplot2::ggplot(df4, ggplot2::aes(x = .data$date, y = .data$logp)) +
      ggplot2::geom_col(fill = "darkblue", color = "darkblue", alpha = 0.75)
  }
  p4 <- p4 +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    panel_theme

  # -------------------- Panel 5: Threshold tiles --------------------
  res_tile <- excode_summary %>%
    dplyr::mutate(
      prob     = as.numeric(1 - excode_summary$posterior0 > posterior_cutoff),
      zscore   = dplyr::if_else(excode_summary$zscore > zscore_cutoff, 2L, 0L),
      pval_bin = dplyr::if_else(excode_summary$pval < p_cutoff, 3L, 0L)
    ) 
  res_tile <- res_tile[,c("date", "prob", "zscore", "pval_bin")] %>%
    tidyr::pivot_longer(cols = 2:4, names_to = "group", values_to = "value") %>%
    dplyr::mutate(
      gg = "Excess\n",
      group = factor(
        dplyr::case_when(
          group == "pval_bin" ~ sprintf("p-value<%s", format(p_cutoff, digits = 3, scientific = FALSE)),
          group == "prob" ~ sprintf("Pr(excess)>%s", format(posterior_cutoff, digits = 3, scientific = FALSE)),
          group == "zscore" ~ sprintf("zscore>%s", format(zscore_cutoff, digits = 3, scientific = FALSE))
        ),
        levels = rev(c(
          sprintf("Pr(excess)>%s", format(posterior_cutoff, digits = 3, scientific = FALSE)),
          sprintf("zscore>%s", format(zscore_cutoff, digits = 3, scientific = FALSE)),
          sprintf("p-value<%s", format(p_cutoff, digits = 3, scientific = FALSE))
        ))
      )
    )

  p5 <- ggplot2::ggplot(res_tile, ggplot2::aes(x = .data$date, y = .data$group, fill = factor(.data$value))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("#00000000", "lightblue", "cornflowerblue", "darkblue")) +
    ggplot2::facet_wrap(~gg, strip.position = "right") +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    panel_theme +
    ggplot2::theme(legend.position = "none")

  # -------------------- Combine --------------------
  p1 + p5 + p2 + p3 + p4 +
    patchwork::plot_layout(heights = c(5, 1, 1.5, 1.5, 1.5), axes = "collect")
}
