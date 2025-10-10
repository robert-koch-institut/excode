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
library(patchwork)
library(tidyr)

## ----load package-------------------------------------------------------------
library(excode)

## ----sts object---------------------------------------------------------------
## Load data and show first 6 rows
data("shadar_df")
kable(as_tibble(shadar_df[1:6, ])) %>%
  kable_styling(font_size = 11)

## ----run_harmonic_model, message=FALSE, warning=FALSE-------------------------
result_shadar_har <- run_excode(
    surv_ts = shadar_df,
    timepoints = 209:295,
    distribution = "Poisson",
    states = 2,
    periodic_model = "Harmonic",
    time_trend = "Linear",
    return_full_model = FALSE, 
    set_baseline_state = TRUE
)

## -----------------------------------------------------------------------------
summary_shadar_har <- summary(result_shadar_har)
kable(as_tibble(summary_shadar_har[82:87, ])) %>%
  kable_styling(font_size = 11)

## ----fig.width=7, fig.height=6, fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 1:** Weekly *Salmonella Hadar* cases with alarm thresholds (upper bounds)."----
#plot_excode_summary(summary_shadar_har, type="bar")
plot_excode_summary(bind_rows(shadar_df %>% 
                                dplyr::filter(!date %in% summary_shadar_har$date & 
                                                lubridate::year(shadar_df$date)>=2003), 
                              summary_shadar_har) %>% 
                                mutate(posterior0=ifelse(is.na(posterior0), 1, posterior0),
                                       zscore=ifelse(is.na(zscore), 0, zscore),
                                        pval=ifelse(is.na(pval), 1, pval))
                    , type="bar")


## ----fig.width=6, fig.height=4------------------------------------------------
data(sarscov2_df)
num_wday <- wday(sarscov2_df$date, week_start = 1)
weekday <- data.frame(day = factor(ifelse(num_wday <= 5, "workday", "weekend"),
  levels = c("workday", "weekend")
))

## ----message=FALSE, warning=FALSE---------------------------------------------

result_sarscov2_custom <- run_excode(
    surv_ts = sarscov2_df,
    covariate_df = weekday,
    states = 2,
    distribution = "Poisson",
    periodic_model = "Custom",
    timepoints = 93:154,
    time_units_back = 8,
    period_length = 7,
    return_full_model = FALSE, 
    set_baseline_state = TRUE
)

summary_sarscov2_custom <- summary(result_sarscov2_custom)



## ----warning=FALSE, message=FALSE, fig.width=7, fig.height=6, fig.align="center", fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 2:** Daily reported SARS-CoV-2 infections in Berlin-Neukölln with alarm thresholds (upper bounds)."----
plot_excode_summary(bind_rows(sarscov2_df %>% dplyr::filter(!date %in%summary_sarscov2_custom$date), summary_sarscov2_custom), type="bar")

## ----data_mort_df_germany, fig.width=6, fig.height=3, fig.align="center", fig.align="center", fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 5:** Weekly all-cause mortaility in Germany."----
data("mort_df_germany")

## ----mort_df_germany, fig.width=6, fig.height=4, fig.align="center"-----------
kable(as_tibble(mort_df_germany[15:20, ])) %>%
  kable_styling(font_size = 11)

## ----multistate_fit, fig.width=7, fig.height=5, fig.align="center"------------

result_mort <- run_excode(
    surv_ts = mort_df_germany,
    covariate_df = weekday,
    states = 3,
    distribution = "NegBinom",
    periodic_model = "Harmonic",
    time_trend = "Spline2",
    timepoints = nrow(mort_df_germany),
    return_full_model = TRUE
)

summary_mort <- summary(result_mort)


## ----multistate_plot, fig.width=7, fig.height=6, fig.align="center", fig.align="center", message=FALSE, warning=FALSE, fig.cap = "**Figure 6:** *Top*: Weekly all-cause mortaility in Germany and expected mortality of the normal and three excess states are shown. Weeks with excess probabilty>0.5 are shown in blue. *Bottom*: Weekly excess probability."----

plot_excode_summary(summary_mort)


## ----echo=FALSE---------------------------------------------------------------
options(old_options)

