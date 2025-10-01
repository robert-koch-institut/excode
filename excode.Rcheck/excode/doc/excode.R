## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_options <- options()
options(digits = 3)
library(knitr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(kableExtra)
library(MASS)
library(gridExtra)

## ----example_mean_har_farr_model, echo=FALSE, fig.cap = "**Figure 1:** Examples for the 'Mean', 'Harmonic' and 'FarringtonNoufaily' model.", fig.align="center", fig.width=6, fig.height=5----
library(excode)
data("snewport_df")
data("shadar_df")
# Use a 'Poisson' distribution in excodeFamily
excode_family_pois <- excodeFamily("Poisson")
# Define the 'Mean' model withouth time trend in excodeFormula
excode_formula_mean <- excodeFormula("Mean",
  timeTrend = FALSE
)
excode_mean_pois <- excodeModel(
  excode_family_pois,
  excode_formula_mean
)
mean_model_pois <- run_excode(snewport_df,
  excode_mean_pois, 410,
  return_full_model = TRUE
)

excode_formula_har <- excodeFormula("Harmonic")
excode_har_pois <- excodeModel(
  excode_family_pois,
  excode_formula_har
)

har_model_pois <- run_excode(shadar_df, excode_har_pois, 285,
  return_full_model = TRUE
)


excode_formula_fn <- excodeFormula("FarringtonNoufaily")
excode_fn_pois <- excodeModel(
  excode_family_pois,
  excode_formula_fn
)

fn_model_pois <- run_excode(shadar_df, excode_fn_pois, 285,
  return_full_model = TRUE
)


fn_df <- data.frame(
  date = fn_model_pois@date,
  observed = fn_model_pois@observed,
  normal = fn_model_pois@emission@mu0,
  excess = fn_model_pois@emission@mu1,
  model = "FarringtonNoufaily"
)
har_df <- data.frame(
  date = har_model_pois@date,
  observed = har_model_pois@observed,
  normal = har_model_pois@emission@mu0,
  excess = har_model_pois@emission@mu1,
  model = "Harmonic"
)
mean_df <- data.frame(
  date = mean_model_pois@date,
  observed = mean_model_pois@observed,
  normal = mean_model_pois@emission@mu0,
  excess = mean_model_pois@emission@mu1,
  model = "Mean"
)

ggplot(bind_rows(mean_df, har_df, fn_df) %>%
  mutate(model = factor(model, levels = c("Mean", "Harmonic", "FarringtonNoufaily")))) +
  geom_col(aes(x = date, y = observed), color = "grey", fill = "lightgrey") +
  geom_line(aes(x = date, y = normal, color = "lightgreen"), linewidth = 1) +
  geom_line(aes(x = date, y = excess, color = "tomato"), linewidth = 1) +
  theme_minimal() +
  ylab("No. of cases") +
  xlab("Week of notification") +
  scale_x_date(date_labels = "%Y", date_breaks = "year") +
  facet_wrap(. ~ model, scales = "free", ncol = 1) +
  scale_colour_manual(
    name = "",
    values = c("lightgreen" = "lightgreen", "tomato" = "tomato"),
    labels = c("normal", "excess")
  ) #+
# theme(legend.position="bottom")

## ----load package-------------------------------------------------------------
library(excode)

## ----sts object---------------------------------------------------------------
## Load data an show first 6 rows
data("shadar_df")
kable(as_tibble(shadar_df[1:6, ])) %>%
  kable_styling(font_size = 11)

## ----plot_shadar, message=FALSE, warning=FALSE, fig.cap = "**Figure 2:** Weekly *Salmonella Newport* cases.", fig.align="center", fig.width=6, fig.height=3----
# Plot the time series
ggplot(shadar_df, aes(x = date, y = observed)) +
  geom_col(color = "grey", fill = "lightgrey") +
  theme_minimal() +
  ylab("No. of cases") +
  xlab("Week of notification") +
  scale_x_date(date_labels = "%Y", date_breaks = "year")

## ----define_poisson_harmonic_model, message=FALSE, warning=FALSE--------------
# Use a 'Poisson' distribution in excodeFamily
excode_family_pois <- excodeFamily("Poisson")
# Define the 'Harmonic' model in excodeFormula
excode_formula_har <- excodeFormula("Harmonic")
# Combine both in the final excodeModel object
excode_har_pois <- excodeModel(
  excode_family_pois,
  excode_formula_har
)

## -----------------------------------------------------------------------------
excode_har_pois

## ----run_harmonic_model-------------------------------------------------------
# Run excode on time points 209-295
result_shadar_har <- run_excode(
  shadar_df,
  excode_har_pois,
  209:295
)

result_shadar_har

## -----------------------------------------------------------------------------
summary_shadar_har <- summary(result_shadar_har)
kable(as_tibble(summary_shadar_har[82:87, ])) %>%
  kable_styling(font_size = 11)

## ----fig.width=6, fig.height=4, fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 3:** Weekly *Salmonella Hadar* cases with alarm thresholds (upper bounds)."----
ggplot(
  shadar_df %>% filter(year(date) >= 2003),
  aes(x = date, y = observed)
) +
  geom_col(color = "grey", fill = "lightgrey") +
  theme_minimal() +
  ylab("No. of cases") +
  xlab("Week of notification") +
  scale_x_date(date_labels = "%Y", date_breaks = "year") +
  geom_step(aes(x = date, y = posterior_ub, color = "upper bound (posterior)"),
    data = summary_shadar_har, linewidth = 1
  ) +
  geom_step(aes(x = date, y = pval_ub, color = "upper bound (pval)"),
    data = summary_shadar_har, linewidth = 1
  ) +
  geom_point(aes(x = date, y = ifelse(observed > pval_ub, observed, NA), color = "Alarm (pval)"),
    data = summary_shadar_har, shape = 3
  ) +
  geom_point(aes(x = date, y = ifelse(observed > posterior_ub, observed, NA), color = "Alarm (posterior)"),
    data = summary_shadar_har, shape = 4
  ) +
  scale_color_manual(values = c("black", "black", "gold", "tomato")) +
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(color = guide_legend(ncol = 2))

## ----fig.width=6, fig.height=4------------------------------------------------
data(sarscov2_df)
num_wday <- wday(sarscov2_df$date, week_start = 1)
weekday <- data.frame(day = factor(ifelse(num_wday <= 5, "workday", "weekend"),
  levels = c("workday", "weekend")
))

## -----------------------------------------------------------------------------
excode_formula_custom <- excodeFormula("Custom", data = weekday, timepoints_per_unit = 7)
excode_custom_pois <- excodeModel(
  excode_family_pois,
  excode_formula_custom
)

## -----------------------------------------------------------------------------
result_sarscov2_custom <- run_excode(sarscov2_df,
  excode_custom_pois,
  93:154,
  time_units_back = 8,
  past_timepoints_not_included = 2
)
summary_sarscov2_custom <- summary(result_sarscov2_custom)

## ----warning=FALSE, message=FALSE, fig.width=6, fig.height=4, fig.align="center", fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 4:** Daily reported SARS-CoV-2 infections in Berlin-Neukölln with alarm thresholds (upper bounds)."----
ggplot(
  sarscov2_df %>%
    filter(date <= max(summary_sarscov2_custom$date)),
  aes(x = date, y = observed)
) +
  geom_col(color = "grey", fill = "lightgrey") +
  theme_minimal() +
  ylab("No. of cases") +
  xlab("Day of reporting") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "month") +
  geom_step(aes(x = date, y = posterior_ub, color = "upper bound (posterior)"),
    data = summary_sarscov2_custom, linewidth = 1
  ) +
  geom_step(aes(x = date, y = pval_ub, color = "upper bound (pval)"),
    data = summary_sarscov2_custom, linewidth = 1
  ) +
  geom_point(aes(x = date, y = ifelse(observed > pval_ub, observed, NA), color = "Alarm (pval)"),
    data = summary_sarscov2_custom, shape = 3
  ) +
  geom_point(aes(x = date, y = ifelse(observed > posterior_ub, observed, NA), color = "Alarm (posterior)"),
    data = summary_sarscov2_custom, shape = 4
  ) +
  scale_color_manual(values = c("black", "black", "gold", "tomato")) +
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(color = guide_legend(ncol = 2))

## ----data_mort_df_germany, fig.width=6, fig.height=3, fig.align="center", fig.align="center", fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 5:** Weekly all-cause mortaility in Germany."----
data("mort_df_germany")

ggplot(mort_df_germany) +
  geom_line(aes(x = date, y = observed), linewidth = 1) +
  theme_minimal() +
  ylab("Number of deaths") +
  xlab("Week of death")

## ----mort_df_germany, fig.width=6, fig.height=4, fig.align="center"-----------
kable(as_tibble(mort_df_germany[1:6, ])) %>%
  kable_styling(font_size = 11)

## ----multistate_bckg_model, fig.width=6, fig.height=4, fig.align="center"-----
# Fit background model
formula_bckg <- "observed ~ sin52 + cos52 + timepoint"
mort_glm_bckg <- glm.nb(formula_bckg, mort_df_germany %>%
  filter(bckg_week))
# Extract prediction for baseline mortality
mort_bckg_mu <- predict(mort_glm_bckg,
  newdata = mort_df_germany,
  type = "response"
)
# store size parameter of ngeative binomial distribution
theta <- mort_glm_bckg$theta

## ----multistate_initial_increase, fig.width=6, fig.height=4, fig.align="center"----
# Relative increase of observed mortality compared to background
p_increase <- mort_df_germany$observed / mort_bckg_mu

# Initial values for multipicative increase of model states
qcut <- c(0.75, 0.85, 0.95)
state_increase <- c(1, quantile(p_increase, prob = qcut))
names(state_increase) <- c("mu0", "mu1", "mu2", "mu3")
state_increase

## ----setup 2, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_options <- options()
options(digits = 6)

## ----initial_mu, fig.width=6, fig.height=4, fig.align="center"----------------
# Initial estimates of expected mortality of all states
initial_mu <- mort_bckg_mu %*% t(state_increase)
kable(as_tibble(initial_mu[1:6, ])) %>%
  kable_styling(font_size = 11)

## ----setup 3, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
old_options <- options()
options(digits = 3)

## ----multistate_model, fig.width=6, fig.height=4, fig.align="center"----------
nStates <- 4
excode_formula_multistate <- excodeFormula("MultiState",
  nStates = nStates,
  intercept = FALSE,
  offset = TRUE
)
excode_family_nb <- excodeFamily("NegBinom", nb_size = theta)
excode_multistate_nb <- excodeModel(excode_family_nb,
  excode_formula_multistate,
  initial_mu = initial_mu
)

## ----multistate_input_data, fig.width=6, fig.height=4, fig.align="center"-----
# Create input data.frame
mort_df_germany_input <- mort_df_germany %>%
  mutate(
    state = ifelse(bckg_week, 0, NA),
    offset = mort_bckg_mu
  ) %>%
  dplyr::select(date, observed, offset, id, state)

## ----multistate_fit, fig.width=7, fig.height=5, fig.align="center"------------
res_excode_multistate_nb <- run_excode(mort_df_germany_input,
  excode_multistate_nb,
  nrow(mort_df_germany_input),
  return_full_model = TRUE,
  time_units_back = 10
)
res_mort_germany <- summary(res_excode_multistate_nb)

## ----multistate_plot, fig.width=7, fig.height=5, fig.align="center", fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 6:** *Top*: Weekly all-cause mortaility in Germany and expected mortality of the normal and three excess states are shown. Weeks with excess probabilty>0.5 are shown in blue. *Bottom*: Weekly excess probability."----
p1 <- ggplot(res_mort_germany %>% mutate(
  obs_excess = ifelse(1 - posterior0 > 0.5, observed, NA),
  group = "Data and model fit"
)) +
  geom_line(aes(x = date, y = observed, color = "Week with excess probabilty<=0.5"),
    linewidth = 1
  ) +
  geom_line(aes(x = date, y = mu0, color = "mu0 (normal)"), linewidth = 0.75) +
  geom_line(aes(x = date, y = mu1, color = "mu1 (low excess)"),
    linewidth = 0.75, alpha = 0.75,
    linetype = "dashed"
  ) +
  geom_line(aes(x = date, y = mu2, color = "mu2 (medium excess)"),
    linewidth = 0.75, alpha = 0.75,
    linetype = "dashed"
  ) +
  geom_line(aes(x = date, y = mu3, color = "mu3 (high excess)"),
    linewidth = 0.75, alpha = 0.75,
    linetype = "dashed"
  ) +
  geom_point(aes(x = date, y = obs_excess, color = "Week with excess probabilty>0.5")) +
  geom_line(aes(x = date, y = obs_excess, color = "Week with excess probabilty>0.5"), linewidth = 1) +
  scale_color_manual(values = c("lightgreen", "gold", "orange", "tomato", "black", "cornflowerblue")) +
  facet_wrap(~group, strip.position = "right") +
  xlab("") +
  ylab("Number of deaths") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5), legend.position = "top",
    legend.title = element_blank(), legend.key.width = unit(1, "cm"),
    plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2)))

p2 <- ggplot(res_mort_germany %>% mutate(excess_prob = 1 - posterior0, group = "Excess probability")) +
  geom_ribbon(aes(x = date, ymin = 0, ymax = excess_prob),
    color = "cornflowerblue", fill = "cornflowerblue",
    alpha = 0.75
  ) +
  geom_hline(yintercept = 0.5, linetype = "dotted", alpha = 0.5) +
  facet_wrap(~group, strip.position = "right") +
  xlab("Week of death") +
  ylab("Excess probability") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"))

layout <- matrix(c(1, 1, 2), ncol = 1)
grid.arrange(grobs = list(p1, p2), nrow = 2, layout_matrix = layout)

## ----create_model_overview, fig.width=7, fig.height=5, fig.align="center", fig.cap = "**Figure 7:** Components of an *excode* Harmonic model using a Poisson distribution."----
har_model_pois <- run_excode(shadar_df, excode_har_pois, 285,
  return_full_model = TRUE
)
plot(har_model_pois)

## ----echo=FALSE---------------------------------------------------------------
options(old_options)

