library(ggplot2)
library(lubridate)
library(surveillance)

sim_df <- NULL

for (n in 1:2) {
  r <- 0.6
  p <- 0.95
  len <- 400
  A <- 1
  alpha <- 3
  beta <- -0.005
  phi <- 0
  frequency <- 1
  K <- log(1.75)

  state <- matrix(data = 0, ncol = 1, nrow = len)
  state[1] <- 0
  transitionMatrix <- matrix(data = c(p, (1 - r), (1 -
    p), r), nrow = 2, ncol = 2)
  if (length(state) > 1) {
    for (i in 2:len) {
      state[i] <- rbinom(1, 1, transitionMatrix[state[i -
        1] + 1, 2])
    }
  }

  state[325:328] <- 1
  observed <- sim.seasonalNoise(
    A, alpha, beta, phi, len,
    frequency, state, K
  )$seasonalBackground


  curr_df <- data.frame(
    date = as_date("2016-01-05") + (1:400) * 7,
    observed = observed,
    id = paste0("sim_", n),
    state = state
  )
  sim_df <- rbind(sim_df, curr_df)
}

ggplot(sim_df, aes(x = date, y = observed)) +
  geom_col(color = "grey") +
  theme_minimal() +
  facet_wrap(. ~ id, nrow = 2)

usethis::use_data(sim_df, internal = TRUE, overwrite = TRUE)
