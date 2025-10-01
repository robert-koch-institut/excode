plot_markov_model <- function(excode_model) {
  transitions <- excode_model@transitions
  init_prob <- excode_model@initial_prob

  par(mar = c(1, 1, 1, 1))
  # Create empty plot
  plot(1:10, type = "n", axes = F, xlab = "", ylab = "")

  # First panel
  rect(1.25, 8.325, 9.975, 10.075)
  text(1.1, 9.9, "Each time point (t) is in", pos = 4, font = 2)
  text(1.1, 9.65, "one of two states (s):", pos = 4, font = 2)
  text(7.5, 9.2, "Normal state (N): ", pos = 2)
  text(7.5, 8.65, "Excess state (E): ", pos = 2)
  symbols(8, 9.2, circles = 0.45, add = TRUE, bg = "springgreen2", fg = "black", inches = FALSE)
  symbols(8, 8.65, circles = 0.45, add = TRUE, bg = "tomato", fg = "black", inches = FALSE)


  # Second panel
  rect(1.25, 6.025, 9.975, 8.275)
  text(1.1, 8.1, "Probability in 1st week", pos = 4, font = 2)
  init_p1 <- round(init_prob[1], 2)
  text(1.1, 7.7, as.expression(bquote("Pr(s"[1] * "=N)=" * .(init_p1))), pos = 4)
  symbols(x = 4.1, y = 7.3, circles = 0.45, add = TRUE, bg = "springgreen2", fg = "black", inches = FALSE)
  arrows(4.55, 7.3, x1 = 6, y1 = 7.3, length = 0.1, lwd = 1, lty = 3)
  init_p2 <- round(init_prob[2], 2)
  text(1.1, 6.7, as.expression(bquote("Pr(s"[1] * "=E)=" * .(init_p2))), pos = 4)
  symbols(x = 4.1, y = 6.3, circles = 0.45, add = TRUE, bg = "tomato", fg = "black", inches = FALSE)
  arrows(4.55, 6.3, x1 = 6, y1 = 6.3, length = 0.1, lwd = 1, lty = 3)


  # Third panel
  rect(1.25, 1.25, 9.975, 5.975)
  text(1.1, 5.8, "State probability in week t", pos = 4, font = 2)
  text(1.1, 5.55, "depends on week t-1", pos = 4, font = 2)

  ypos <- 4.6
  ydiff <- 1
  trans_cols <- list(
    c("springgreen2", "springgreen2"),
    c("springgreen2", "tomato"),
    c("tomato", "springgreen2"),
    c("tomato", "tomato")
  )
  transition_names <- list(c("N", "N"), c("E", "N"), c("N", "E"), c("E", "E"))
  for (k in 1:length(trans_cols)) {
    arrows(2.2, ypos, x1 = 3.65, y1 = ypos, length = 0.1, lwd = 1, lty = 3)
    symbols(
      x = 4.1, y = ypos, circles = 0.45, add = TRUE,
      bg = trans_cols[[k]][1], fg = "black", inches = FALSE
    )
    arrows(4.55, ypos, x1 = 6, y1 = ypos, length = 0.1, lwd = 1, lty = 1)
    symbols(
      x = 6.45, y = ypos, circles = 0.45, add = TRUE,
      bg = trans_cols[[k]][2], fg = "black", inches = FALSE
    )
    arrows(6.9, ypos, x1 = 8.45, y1 = ypos, length = 0.1, lwd = 1, lty = 3)
    curr_trans <- round(c(transitions[1, ], transitions[2, ])[k], 2)
    text(1.1, ypos + 0.4, as.expression(bquote("Pr(s"[t] * "=" * .(transition_names[[k]][1]) * "|s"["t-1"] * "=" * .(transition_names[[k]][2]) * ")=" * .(curr_trans))), pos = 4)
    # text(1.1, ypos+0.4, expression(phantom('Pr(s'[t]*'=')*"N"*phantom('|s'[t-1]*'=')*phantom("N")*phantom(')=0.96')), pos=4, col=trans_cols[[k]][1], font=2)
    ypos <- ypos - ydiff
  }
}



plot_glm <- function(excode_model,
                     set_mar = TRUE,
                     cols = c("lightgrey", "grey", "blue"),
                     state_cols = c("springgreen2", "tomato")) {
  y <- excode_model@observed
  pred <- matrix(c(
    excode_model@emission@mu0,
    excode_model@emission@mu1
  ), ncol = 2)
  surv_ts <- sts(
    observed = y, epoch = excode_model@date,
    frequency = excode_model@emission@excode_formula@timepoints_per_unit
  )

  # plot glm
  if (set_mar) {
    par(mar = c(2, 2, 3, 0.5) + 0.1)
  }
  yMax <- max(c(pred, y))
  yLim <- c(-yMax * 0.05, yMax * 1.15)
  legend.offset <- 0
  legend.offset_y <- 0
  plot(surv_ts,
    main = "", col = cols, legend = F,
    xlab = "Time (weeks)",
    xaxis.tickFreq = list("%m" = atChange, "%Y" = atChange),
    xaxis.labelFreq = list("%Y" = atMedian), xaxis.labelFormat = "%Y",
    outbreak.symbol = list(pch = 3, col = "red", cex = 1, lwd = 1), ylim = yLim
  )


  lines(pred[, 1], col = state_cols[1], lwd = 2)
  lines(pred[, 2], col = state_cols[2], lwd = 2)


  legend("top", c("Expected cases:", "", "Normal", "Excess"),
    ncol = 2,
    col = c("#00000000", "#00000000", state_cols[1], state_cols[2]), lwd = 2, bty = "n"
  )
}


plot_posterior_result <- function(excode_model) {
  cols <- c("springgreen2", "tomato")
  par(mar = c(3, 1, 2, 2) + 0.1)

  # plot hidden states
  plot(1:10, xaxt = "n", yaxt = "n", ylab = "", xlab = "", type = "n")
  text(5.5, 8.5, "Alarm based on Posterior", pos = 3, font = 2)
  circle_col <- ifelse(excode_model@posterior[length(excode_model@posterior)] > 0.5, cols[2], cols[1])
  symbols(x = 9 + 0.5, y = 5, circles = 0.3, add = TRUE, bg = circle_col, fg = "black", inches = FALSE)
  arrows(7.25 + 0.5, 5, x1 = 8.7 + 0.5, y1 = 5, length = 0.1, lwd = 1)

  circle_col <- ifelse(excode_model@posterior[length(excode_model@posterior) - 1] > 0.5, cols[2], cols[1])
  symbols(x = 6.95 + 0.5, y = 5, circles = 0.3, add = TRUE, bg = circle_col, fg = "black", inches = FALSE)
  arrows(5.2 + 0.5, 5, x1 = 6.65 + 0.5, y1 = 5, length = 0.1, lwd = 1)

  circle_col <- ifelse(excode_model@posterior[length(excode_model@posterior) - 2] > 0.5, cols[2], cols[1])
  symbols(x = 4.9 + 0.5, y = 5, circles = 0.3, add = TRUE, bg = circle_col, fg = "black", inches = FALSE)
  arrows(3.15 + 0.5, 5, x1 = 4.6 + 0.5, y1 = 5, length = 0.1, lwd = 1, lty = 3)

  # plot observed data
  arrows(9 + 0.5, 5 + 0.3 + 0.15, x1 = 9 + 0.5, y1 = 5 + 1.75, length = 0.1, lwd = 1)
  text(9 + 0.5, 5 + 1.75, excode_model@observed[length(excode_model@observed)], pos = 3)
  arrows(6.95 + 0.5, 5 + 0.3 + 0.15, x1 = 6.95 + 0.5, y1 = 5 + 1.75, length = 0.1, lwd = 1)
  text(6.95 + 0.5, 5 + 1.75, excode_model@observed[length(excode_model@observed) - 1], pos = 3)
  arrows(4.9 + 0.5, 5 + 0.3 + 0.15, x1 = 4.9 + 0.5, y1 = 5 + 1.75, length = 0.1, lwd = 1)
  text(4.9 + 0.5, 5 + 1.75, excode_model@observed[length(excode_model@observed) - 2], pos = 3)

  # plot description
  text(0.75, 5 + 1.75 + 0.75, "No. of cases", pos = 4, offset = 0)
  text(0.75, 5 + 0.25, "States", pos = 4, offset = 0)
  text(0.75, 3.5 + 0.25, "Time point", pos = 4, offset = 0)
  text(9 + 0.5, 3.5, excode_model@timepoint[length(excode_model@timepoint)], pos = 3, offset = 0)
  text(6.95 + 0.5, 3.5, excode_model@timepoint[length(excode_model@timepoint) - 1], pos = 3, offset = 0)
  text(4.9 + 0.5, 3.5, excode_model@timepoint[length(excode_model@timepoint) - 2], pos = 3, offset = 0)

  obs_data <- excode_model@observed[(length(excode_model@observed) - 2):length(excode_model@observed)]
  # plot posterior descriotion (O)
  obs_text <- paste0(paste(rev(obs_data), collapse = ","), ", ...")
  obs_text <- paste0("=E|", obs_text, "):")
  add_text <- as.expression(bquote("Pr(s"["t"] * .(obs_text)))
  text(0.75, 1 + 0.5, add_text, pos = 4, offset = 0)
  # plot posterior descriotion (E)
  obs_text <- paste0(paste(rev(obs_data), collapse = ","), ", ...")
  obs_text <- paste0("=N|", obs_text, "):")
  add_text <- as.expression(bquote("Pr(s"["t"] * .(obs_text)))
  text(0.75, 2 + 0.5, add_text, pos = 4, offset = 0)


  num <- sprintf("%.4f", round(excode_model@posterior[length(excode_model@posterior)], 4))
  hh <- as.expression(bquote(bold(.(num))))
  text(9 + 0.95, 1 + 0.5, hh, pos = 2, offset = 0, font = 2)
  num <- sprintf("%.4f", round(1 - round(excode_model@posterior[length(excode_model@posterior)], 4), 4))
  hh <- as.expression(bquote(.(num)))
  text(9 + 0.95, 2 + 0.5, hh, pos = 2, offset = 0, font = 2)

  num <- sprintf("%.4f", round(excode_model@posterior[length(excode_model@posterior) - 1], 4))
  hh <- as.expression(bquote(.(num)))
  text(6.95 + 0.95, 1 + 0.5, hh, pos = 2, offset = 0, font = 2)
  num <- sprintf("%.4f", round(1 - round(excode_model@posterior[length(excode_model@posterior) - 1], 4), 4))
  hh <- as.expression(bquote(.(num)))
  text(6.95 + 0.95, 2 + 0.5, hh, pos = 2, offset = 0, font = 2)

  num <- sprintf("%.4f", round(excode_model@posterior[length(excode_model@posterior) - 2], 4))
  hh <- as.expression(bquote(.(num)))
  text(4.9 + 0.95, 1 + 0.5, hh, pos = 2, offset = 0, font = 2)
  num <- sprintf("%.4f", round(1 - round(excode_model@posterior[length(excode_model@posterior) - 2], 4), 4))
  hh <- as.expression(bquote(.(num)))
  text(4.9 + 0.95, 2 + 0.5, hh, pos = 2, offset = 0, font = 2)
}


plot_pval_result <- function(excode_model) {
  par(mar = c(3, 4, 2, 2) + 0.1)
  upper <- summary(excode_model)$pval_ub
  upper <- upper[length(upper)]
  xpos <- 0:(round(max(excode_model@observed) * 1.2))

  index <- which(excode_model@timepoint_fit == excode_model@timepoint)

  curr_data <- create_emission_prob_input(
    excode_model@emission@distribution,
    xpos,
    rep(excode_model@emission@mu0[index], length(xpos)),
    rep(excode_model@emission@mu1[index], length(xpos)),
    index
  )
  curr_data <- curr_data[curr_data$state == 0, ]
  prob <- calcEmissionProb(
    excode_model@emission@distribution,
    curr_data
  )[, 1]


  timepoint <- max(excode_model@timepoint)
  plot(xpos, prob,
    type = "n", ylab = as.expression(bquote("Pr(#cases | s"[.(timepoint)] * "=E)")),
    xlab = "No. of cases", ylim = c(0, max(prob * 1.35))
  )
  text(max(xpos / 2), max(prob * 1.35), "Alarm based\non p-value", pos = 1, font = 2)

  for (i in 1:(length(xpos) - 1)) {
    rect(xpos[i], 0, xpos[i + 1], prob[i], col = "springgreen2")
  }
  abline(v = upper, lty = 2)
  pval <- summary(excode_model)$pval
  pval <- pval[length(pval)]
  axis(3, at = upper, paste0("p-value=", sprintf("%.6f", pval)), font = 2)
}



#' Plot sts object with excode results
#'
#' @title Plot sts object with excode results
#'
#' @param excode_model An object of class \code{\linkS4class{excodeModel}} which specifies the model parameters and structure
#'
#' @examples
#'
#' # TODO
#'
#' @export
plot_model <- function(excode_model) {
  layout_mat <- matrix(c(
    1, 1, 2, 2, 2, 2, 2,
    1, 1, 3, 3, 3, 4, 4
  ), byrow = T, ncol = 7)
  layout(layout_mat,
    widths = c(1, 1, 1),
    heights = c(1, 1)
  )
  plot_markov_model(excode_model)

  plot_glm(excode_model)

  plot_posterior_result(excode_model)

  plot_pval_result(excode_model)

  dist <- excode_model@emission@distribution@name
  mod <- excode_model@emission@excode_formula@name
  title(paste0(
    "Model fitted at time point ", max(excode_model@timepoint), " with model ",
    dist, "-", mod
  ), outer = T, line = -1)
}
