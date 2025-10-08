#' Plot excode summary
#'
#' Creates a multi-panel diagnostic plot summarizing excess mortality model results,
#' including observed deaths, model fits, posterior excess probability, Anscombe residuals,
#' and p-values. Allows user-defined cutoff thresholds for each diagnostic.
#'
#' @param excode_summary A data frame containing model summary output with columns:
#'   \code{date}, \code{observed}, \code{posterior0}, \code{anscombe_residual},
#'   \code{pval}, and one or more fitted mean columns named \code{mu0}, \code{mu1}, etc.
#' @param states Integer. Number of model states (i.e., fitted mean columns \code{mu0}, \code{mu1}, ...).
#'   If missing, this will be inferred automatically from the column names.
#' @param posterior_cutoff Numeric (default = 0.5). Threshold for classifying excess probability
#'   (\eqn{1 - posterior0 > posterior_cutoff}).
#' @param anscombe_cutoff Numeric (default = 2). Threshold for Anscombe residuals.
#' @param p_cutoff Numeric (default = 0.05). Significance threshold for p-values.
#'
#' @return A combined ggplot2 object with multiple panels visualizing:
#' \itemize{
#'   \item Observed deaths and model fits
#'   \item Posterior excess probability
#'   \item Anscombe residuals
#'   \item -log10(p-values)
#'   \item Threshold indicator tiles
#' }
#'
#' @details
#' The function visualizes multiple diagnostics side-by-side using the \pkg{patchwork} layout system.
#' Dotted horizontal lines indicate the chosen thresholds for posterior probability, residuals,
#' and p-values. The tile panel summarizes binary exceedances across the thresholds.
#'
#' @examples
#' \dontrun{
#' plot_excode_summary(
#'   excode_summary,
#'   posterior_cutoff = 0.6,
#'   anscombe_cutoff  = 2.5,
#'   p_cutoff         = 0.01
#' )
#' }
#'
#' @export
#'
plot_excode_summary <- function(
    excode_summary,
    states,
    posterior_cutoff = 0.5,
    anscombe_cutoff  = 2,
    p_cutoff         = 0.05
) {
  # infer number of state curves + mu names if not supplied
  mu_names <- grep("^mu\\d+$", names(excode_summary), value = TRUE)
  if (missing(states)) {
    states <- length(mu_names)
  } else {
    mu_names <- paste0("mu", 0:(states - 1))
    mu_names <- intersect(mu_names, names(excode_summary))
  }
  
  # ---------- Panel 1: observed + state fits (with legend on top) ----------
  # Build long data for observed + all mu*
  df_obs <- excode_summary |>
    dplyr::select(date, observed) |>
    dplyr::rename(value = observed) |>
    dplyr::mutate(series = "Observed")
  
  df_mus <- excode_summary |>
    dplyr::select(date, dplyr::all_of(mu_names)) |>
    tidyr::pivot_longer(-date, names_to = "series", values_to = "value")
  
  df_p1 <- dplyr::bind_rows(df_obs, df_mus) |>
    dplyr::mutate(group = "Number of deaths\nand model fit")
  
  # Color palette: keep original style (greens for mu0, gold→tomato for others, black for observed)
  mu_cols <- c(
    grDevices::colorRampPalette(c("lightgreen", "forestgreen"))(5)[2],
    if (length(mu_names) > 1) grDevices::colorRampPalette(c("gold", "tomato"))(length(mu_names) - 1) else character(0)
  )
  names(mu_cols) <- mu_names
  legend_values <- c("Observed" = "black", mu_cols)
  
  p1 <- ggplot2::ggplot(df_p1, ggplot2::aes(x = date, y = value, color = series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = legend_values, breaks = names(legend_values)) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("Week of death") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.key.width = grid::unit(1, "cm"),
      plot.margin = grid::unit(c(0, 0.5, 0, 0), "cm")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1, byrow = TRUE))
  
  # ---------- Panel 2: posterior excess probability ----------
  p2 <- ggplot2::ggplot(excode_summary |>
                          dplyr::mutate(excess_prob = 1 - posterior0, group = "Excess\nProbability")) +
    ggplot2::geom_ribbon(ggplot2::aes(x = date, ymin = 0, ymax = excess_prob),
                         color = "lightblue", fill = "lightblue", alpha = 0.75) +
    ggplot2::geom_hline(yintercept = posterior_cutoff, linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("Week of death") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(breaks = c(0, posterior_cutoff, 1)) +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0.5, 0, 0), "cm"))
  
  # ---------- Panel 3: Anscombe residuals ----------
  p3 <- ggplot2::ggplot(excode_summary |>
                          dplyr::mutate(group = "Anscombe\nresidual")) +
    ggplot2::geom_ribbon(ggplot2::aes(x = date, ymin = 0, ymax = anscombe_residual),
                         color = "cornflowerblue", fill = "cornflowerblue", alpha = 0.75) +
    ggplot2::geom_hline(yintercept = anscombe_cutoff, linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("Week of death") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0.5, 0, 0), "cm"))
  
  # ---------- Panel 4: -log10(p) ----------
  p4 <- ggplot2::ggplot(excode_summary |>
                          dplyr::mutate(group = "-log10(pval)\n")) +
    ggplot2::geom_ribbon(ggplot2::aes(x = date, ymin = 0, ymax = -log10(pval)),
                         color = "darkblue", fill = "darkblue", alpha = 0.75) +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), linetype = "dotted", alpha = 0.5) +
    ggplot2::facet_wrap(~group, strip.position = "right") +
    ggplot2::xlab("Week of death") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0.5, 0, 0), "cm"))
  
  # ---------- Panel 5: binary threshold tiles ----------
  res_tile <- excode_summary |>
    dplyr::mutate(
      prob     = as.numeric(1 - posterior0 > posterior_cutoff),
      anscombe = ifelse(anscombe_residual > anscombe_cutoff, 2, 0),
      pval_bin = ifelse(pval < p_cutoff, 3, 0)
    ) |>
    dplyr::select(date, prob, anscombe, pval_bin) |>
    tidyr::pivot_longer(cols = 2:4, names_to = "group", values_to = "value") |>
    dplyr::mutate(
      gg = "Excess\n",
      group = factor(
        dplyr::case_when(
          group == "pval_bin" ~ sprintf("p-value<%s", format(p_cutoff, digits = 3, scientific = FALSE)),
          group == "prob"     ~ sprintf("Pr(excess)>%s", format(posterior_cutoff, digits = 3, scientific = FALSE)),
          group == "anscombe" ~ sprintf("AR>%s", format(anscombe_cutoff, digits = 3, scientific = FALSE))
        ),
        levels = rev(c(
          sprintf("Pr(excess)>%s", format(posterior_cutoff, digits = 3, scientific = FALSE)),
          sprintf("AR>%s",        format(anscombe_cutoff,  digits = 3, scientific = FALSE)),
          sprintf("p-value<%s",   format(p_cutoff,         digits = 3, scientific = FALSE))
        ))
      )
    )
  
  p5 <- ggplot2::ggplot(res_tile, ggplot2::aes(x = date, y = group, fill = factor(value))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("white", "lightblue", "cornflowerblue", "darkblue")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Week of death") + ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~gg, strip.position = "right") +
    ggplot2::theme(plot.margin = grid::unit(c(0, 0.5, 0, 0), "cm"),
                   legend.position = "none")
  
  # ---------- combine panels ----------
  p1 + p5 + p2 + p3 + p4 +
    patchwork::plot_layout(heights = c(5, 1, 1.5, 1.5, 1.5), axes = "collect")
}
