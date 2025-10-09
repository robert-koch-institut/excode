#' Plot excode summary
#'
#' Creates a multi-panel diagnostic plot summarizing excess mortality model results,
#' including observed cases, model fits, posterior excess probability, Anscombe residuals,
#' and p-values. Allows user-defined cutoff thresholds for each diagnostic.
#'
#' @param excode_summary A data frame containing model summary output with columns:
#'   \code{date}, \code{observed}, \code{posterior0}, \code{zscore},
#'   \code{pval}, and one or more fitted mean columns named \code{mu0}, \code{mu1}, etc.
#' @param states Integer. Number of model states (i.e., fitted mean columns \code{mu0}, \code{mu1}, ...).
#'   If missing, this will be inferred automatically from the column names.
#' @param posterior_cutoff Numeric (default = 0.5). Threshold for classifying excess probability
#'   (\eqn{1 - posterior0 > posterior_cutoff}).
#' @param zscore_cutoff Numeric (default = 2.33). Threshold for Anscombe residuals.
#' @param p_cutoff Numeric (default = 0.01). Significance threshold for p-values.
#' @param type Character. Either "line" (default) or "bar".
#'   - Top panel: mu* always lines/steps; observed is line ("line") or bars ("bar").
#'   - Diagnostics: ribbons ("line") or bars ("bar").
#'
#' @return A combined ggplot2 object.
#' @export
plot_excode_summary <- function(
    excode_summary,
    states,
    posterior_cutoff = 0.5,
    zscore_cutoff  = 2.33,
    p_cutoff         = 0.01,
    type             = c("line", "bar")
) {
  type <- match.arg(type)
  
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
  df_obs <- excode_summary |>
    dplyr::select(date, observed) |>
    dplyr::rename(value = observed) |>
    dplyr::mutate(series = "Observed")
  
  df_mus <- excode_summary |>
    dplyr::select(date, dplyr::all_of(mu_names)) |>
    tidyr::pivot_longer(-date, names_to = "series", values_to = "value")
  
  df_p1 <- dplyr::bind_rows(df_obs, df_mus) |>
    dplyr::mutate(group = "Number of cases\nand model fit")
  
  # Color palette: black observed; greens for mu0, gold→tomato for remaining mu
  mu_cols <- c(
    grDevices::colorRampPalette(c("lightgreen", "forestgreen"))(5)[2],
    if (length(mu_names) > 1) grDevices::colorRampPalette(c("gold", "tomato"))(length(mu_names) - 1) else character(0)
  )
  names(mu_cols) <- mu_names
  legend_values <- c("Observed" = "black", mu_cols)
  
  p1 <- ggplot2::ggplot(df_p1, ggplot2::aes(x = date, y = value, color = series))
  if (type == "bar") {
    # Observed as bars; mu* as steps
    p1 <- p1 +
      ggplot2::geom_col(
        data = subset(df_p1, series == "Observed"),
        fill = "lightgrey", color = "grey", alpha = 0.9
      ) +
      ggplot2::geom_step(
        data = subset(df_p1, series != "Observed"),
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
  df2 <- excode_summary |>
    dplyr::mutate(excess_prob = 1 - posterior0, group = "Excess\nProbability")
  
  if (type == "line") {
    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = excess_prob),
                           fill = "lightblue", color = "lightblue", alpha = 0.75)
  } else {
    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = date, y = excess_prob)) +
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
  df3 <- excode_summary |>
    dplyr::mutate(
      ymin_zscore = pmin(0, zscore),
      ymax_zscore = pmax(0, zscore),
      group   = "zscore\n"
    )
  
  if (type == "line") {
    p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin_zscore, ymax = ymax_zscore),
                           fill = "cornflowerblue", color = "cornflowerblue", alpha = 0.75)
  } else {
    p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = date, y = zscore)) +
      ggplot2::geom_col(fill = "cornflowerblue", color = "cornflowerblue", alpha = 0.75)
  }
  p3 <- p3 +
    ggplot2::geom_hline(yintercept = zscore_cutoff, linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    panel_theme
  
  # -------------------- Panel 4: -log10(pval) --------------------
  df4 <- excode_summary |>
    dplyr::mutate(logp = -log10(pval), group = "-log10(pval)\n")
  
  if (type == "line") {
    p4 <- ggplot2::ggplot(df4, ggplot2::aes(x = date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = logp),
                           fill = "darkblue", color = "darkblue", alpha = 0.75)
  } else {
    p4 <- ggplot2::ggplot(df4, ggplot2::aes(x = date, y = logp)) +
      ggplot2::geom_col(fill = "darkblue", color = "darkblue", alpha = 0.75)
  }
  p4 <- p4 +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    panel_theme
  
  # -------------------- Panel 5: Threshold tiles --------------------
  res_tile <- excode_summary |>
    dplyr::mutate(
      prob     = as.numeric(1 - posterior0 > posterior_cutoff),
      zscore = ifelse(zscore > zscore_cutoff, 2, 0),
      pval_bin = ifelse(pval < p_cutoff, 3, 0)
    ) |>
    dplyr::select(date, prob, zscore, pval_bin) |>
    tidyr::pivot_longer(cols = 2:4, names_to = "group", values_to = "value") |>
    dplyr::mutate(
      gg = "Excess\n",
      group = factor(
        dplyr::case_when(
          group == "pval_bin" ~ sprintf("p-value<%s", format(p_cutoff, digits = 3, scientific = FALSE)),
          group == "prob"     ~ sprintf("Pr(excess)>%s", format(posterior_cutoff, digits = 3, scientific = FALSE)),
          group == "zscore" ~ sprintf("zscore>%s",        format(zscore_cutoff,  digits = 3, scientific = FALSE))
        ),
        levels = rev(c(
          sprintf("Pr(excess)>%s", format(posterior_cutoff, digits = 3, scientific = FALSE)),
          sprintf("zscore>%s",        format(zscore_cutoff,  digits = 3, scientific = FALSE)),
          sprintf("p-value<%s",   format(p_cutoff,         digits = 3, scientific = FALSE))
        ))
      )
    )
  
  p5 <- ggplot2::ggplot(res_tile, ggplot2::aes(x = date, y = group, fill = factor(value))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("#00000000", "lightblue", "cornflowerblue", "darkblue")) +
    ggplot2::facet_wrap(~gg, strip.position = "right") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    panel_theme +
    ggplot2::theme(legend.position = "none")
  
  # -------------------- Combine --------------------
  p1 + p5 + p2 + p3 + p4 +
    patchwork::plot_layout(heights = c(5, 1, 1.5, 1.5, 1.5), axes = "collect")
}
