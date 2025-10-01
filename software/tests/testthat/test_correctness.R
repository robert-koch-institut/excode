library(excode)
library(testthat)
library(dplyr)
suppressWarnings(library(surveillance))

# retrieving the internal datasets from sysdata.rda
sim_df <- get0("sim_df", envir = asNamespace("excode"))
test_results <- get0("test_results", envir = asNamespace("excode"))

multiple_ts <- sim_df


## Create models for testing
## Negative Binomial
# Harmonic
excode_har_nb <- excodeModel(
  excodeFamily("NegBinom"),
  excodeFormula("Harmonic")
)
# Farrington-Noufaily
excode_fn_nb <- excodeModel(
  excodeFamily("NegBinom"),
  excodeFormula("FarringtonNoufaily", w = 4)
)
## Poisson
# Harmonic
excode_har_pois <- excodeModel(
  excodeFamily("Poisson"),
  excodeFormula("Harmonic")
)
# Farrington-Noufaily
excode_fn_pois <- excodeModel(
  excodeFamily("Poisson"),
  excodeFormula("FarringtonNoufaily", w = 4)
)



data_har <- excode:::prepareData(single_ts, excode_har_nb, 325,
  id = "sim_1", time_units_back = 5,
  past_weeks_not_included_state = 26,
  past_weeks_not_included_init = 26
)

data_har_custom <- bind_rows(
  data_har %>% mutate(id = "sim_1"),
  data_har %>% mutate(id = "sim_2")
)

var_sel <- c("timepoint", "sin1", "cos1", "id")
# NegBinom
excode_custom_har_nb_multi <- excodeModel(
  excodeFamily("NegBinom"),
  excodeFormula("Custom",
    data = data_har_custom[, var_sel]
  )
)
var_sel <- c("timepoint", "sin1", "cos1")
excode_custom_har_nb_single <- excodeModel(
  excodeFamily("NegBinom"),
  excodeFormula("Custom",
    data = data_har[, var_sel]
  )
)

# Poisson
excode_custom_har_pois_multi <- excodeModel(
  excodeFamily("Poisson"),
  excodeFormula("Custom",
    data = data_har_custom[, var_sel]
  )
)

var_sel <- c("timepoint", "sin1", "cos1")
excode_custom_har_pois_single <- excodeModel(
  excodeFamily("Poisson"),
  excodeFormula("Custom",
    data = data_har[, var_sel]
  )
)


## Test models
## Negative Binomial
# Harmonic
res_har_single_nb <- list(
  summary(run_excode(surv_ts=single_ts[,1:3], excode_model=excode_har_nb, timepoints=325)),
  summary(run_excode(surv_ts=single_ts, excode_model=excode_har_nb, timepoints=325))
)
# Farrington Noufaily
res_fn_single_nb <- list(
  summary(run_excode(surv_ts=single_ts[,1:3], excode_model=excode_fn_nb, timepoints=325)),
  summary(run_excode(surv_ts=single_ts, excode_model=excode_fn_nb, timepoints=325))
)
# Custom
res_custom_single_nb <- list(
  summary(run_excode(surv_ts=single_ts[data_har$rtime, 1:3], excode_model=excode_custom_har_nb_single, length(data_har$rtime))),
  summary(run_excode(surv_ts=single_ts[data_har$rtime, ], excode_model=excode_custom_har_nb_single, length(data_har$rtime)))
)




## Poisson
# Harmonic
res_har_single_pois <- list(
  summary(run_excode(surv_ts=single_ts[,1:3], excode_model=excode_har_pois, timepoints=325)),
  summary(run_excode(surv_ts=single_ts, excode_model=excode_har_pois, timepoints=325))
)
# Farrington Noufaily
res_fn_single_pois <- list(
  summary(run_excode(surv_ts=single_ts[,1:3], excode_model=excode_fn_pois, timepoints=325)),
  summary(run_excode(surv_ts=single_ts, excode_model=excode_fn_pois, timepoints=325))
)
# Custom
res_custom_single_pois <- list(
  summary(run_excode(surv_ts=single_ts[data_har$rtime,1:3], excode_model=excode_custom_har_pois_single, timepoints=length(data_har$rtime))),
  summary(run_excode(surv_ts=single_ts[data_har$rtime, ], excode_model=excode_custom_har_pois_single, timepoints=length(data_har$rtime)))
)





if (FALSE) {
  test_results <- list(
    res_har_single_nb = res_har_single_nb,
    res_har_mult_nb = res_har_mult_nb,
    res_fn_single_nb = res_fn_single_nb,
    res_fn_mult_nb = res_fn_mult_nb,
    res_custom_single_nb = res_custom_single_nb,
    res_custom_mult_nb = res_custom_mult_nb,
    res_har_single_pois = res_har_single_pois,
    res_har_mult_pois = res_har_mult_pois,
    res_fn_single_pois = res_fn_single_pois,
    res_fn_mult_pois = res_fn_mult_pois,
    res_custom_single_pois = res_custom_single_pois,
    res_custom_mult_pois = res_custom_mult_pois
  )
  test_results <- do.call("rbind", lapply(names(test_results), function(x) {
    out <- do.call("rbind", test_results[[x]])[, test_var_pois]
    out$name <- x
    out
  }))
  save(test_results, file = "test_results.rda")
}

test_var_pois <- c("posterior", "pval", "mu0", "mu1", "BIC", "AIC")

all_results <- list(
  res_har_single_nb = res_har_single_nb,
  res_fn_single_nb = res_fn_single_nb,
  res_custom_single_nb = res_custom_single_nb,
  res_har_single_pois = res_har_single_pois,
  res_fn_single_pois = res_fn_single_pois,
  res_custom_single_pois = res_custom_single_pois
)




all_results <- do.call("rbind", lapply(names(all_results), function(x) {
  out <- do.call("rbind", all_results[[x]])[, test_var_pois]
  out$name <- x
  out
}))

#take <- c(1:2, 4:7, 10:11, 13:16, 19:20, 22:25, 28:29, 31:34, 37:38, 40:43, 46:47, 49:52)
#take <- c(1:2, 10:11, 19:20, 28:29, 37:38, 46:47)
#rownames(all_results) <- take
test_that("Consistensy of results", {
  expect_equal(all_results, test_results)
})

