# Summary of an excodeModel

Summarize a fitted `excodeModel` at time points where excess-count
detection was performed, including expectations, posterior
probabilities, p-values, and corresponding thresholds per the model
specification.

## Usage

``` r
# S4 method for class 'excodeModel'
summary(
  object,
  pars = c("posterior", "pval", "zscore", "date", "timepoint", "observed", "emission",
    "BIC", "AIC"),
  prob_threshold = 0.5,
  pval_threshold = 0.05,
  anscombe_threshold = 2,
  maxiter = 1000
)
```

## Arguments

- object:

  An `excodeModel` to summarize.

- pars:

  Character vector of fields to include (e.g., `"posterior"`, `"pval"`,
  `"zscore"`, `"date"`, `"timepoint"`, `"observed"`, `"emission"`,
  `"BIC"`, `"AIC"`); defaults to all.

- prob_threshold:

  Numeric posterior probability threshold used to compute the
  posterior-based upper bound (`posterior_ub`) and flag excess counts
  (\>= threshold); default 0.5.

- pval_threshold:

  Numeric p-value threshold used to compute the p-value-based upper
  bound (`pval_ub`) and flag excess counts (\<= threshold); default
  0.01.

- anscombe_threshold:

  Numeric threshold on the Anscombe z-score used for signal assessment;
  default 2.

- maxiter:

  Integer maximum iterations when estimating the posterior alarm
  threshold; default 1000.

## Value

A `data.frame` summarizing selected components of the `excodeModel`
(expected values, posterior, p-value, and fit metrics such as AIC/BIC)
according to `pars`.

## See also

[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel-class.md),
[`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily-class.md),
[`excodeFormula`](https://robert-koch-institut.github.io/excode/reference/excodeFormula-class.md)

## Examples

``` r
data(shadar_df)
res_har_pois <- run_excode(surv_ts = shadar_df, timepoints = 295, distribution = "Poisson", 
states = 2, periodic_model = "Harmonic", time_trend = "Linear", set_baseline_state = TRUE)
summary(res_har_pois)
#>   posterior0 posterior1      pval   zscore       date timepoint observed
#> 1  0.8669206  0.1330794 0.1883469 1.141199 2006-08-28       295        6
#>        mu0      mu1      BIC      AIC
#> 1 3.826315 12.43947 1149.244 1066.034
```
