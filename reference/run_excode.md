# Detect excess counts in epidemiological time series

Fits an
[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel-class.md)
to surveillance data and returns a fitted model with posterior state
probabilities, fitted means, and diagnostics for excess-count detection;
supports a single time series provided as a `data.frame`.

## Usage

``` r
run_excode(
  surv_ts,
  timepoints = NULL,
  time_units_back = 5,
  distribution,
  states,
  time_trend = "None",
  periodic_model = "Harmonic",
  period_length = 52,
  intercept = TRUE,
  covariate_df = NULL,
  weights_threshold = 2.58,
  weights_threshold_baseline = 1,
  set_baseline_state = FALSE,
  maxIter = 100,
  verbose = FALSE,
  return_full_model = FALSE
)
```

## Arguments

- surv_ts:

  `data.frame` surveillance time series as described in *Details*.

- timepoints:

  Integer or integer vector of time points at which to perform
  detection; must be provided unless embedded in a `MultiState` formula.

- time_units_back:

  Integer number of past time *units* used for fitting (interpreted
  relative to `timepoints_per_unit` in the model’s `excodeFormula`);
  default `5`.

- distribution:

  Character emission family passed to
  [`excodeFamily()`](https://robert-koch-institut.github.io/excode/reference/excodeFamily.md),
  e.g., `"Poisson"` or `"NegBinom"`.

- states:

  Integer (\>= 2) number of latent states for initialization/fitting.

- time_trend:

  Character time-trend specification forwarded to
  [`excodeFormula()`](https://robert-koch-institut.github.io/excode/reference/excodeFormula.md),
  one of `c("Linear","Spline1","Spline2","Spline3","Spline4","None")`;
  default `"None"`.

- periodic_model:

  Character background model name forwarded to
  [`excodeFormula()`](https://robert-koch-institut.github.io/excode/reference/excodeFormula.md),
  one of `c("Mean","FarringtonNoufaily","Harmonic","Custom")`; default
  `"Harmonic"`.

- period_length:

  Integer number of time points per unit season used as
  `timepoints_per_unit` in
  [`excodeFormula()`](https://robert-koch-institut.github.io/excode/reference/excodeFormula.md)
  (e.g., `52` weekly per year, `365` daily per year); default `52`.

- intercept:

  Logical; include an intercept in the constructed formulas; default
  `TRUE`.

- covariate_df:

  `data.frame` or `NULL`; optional covariates used when constructing
  [`excodeFormula()`](https://robert-koch-institut.github.io/excode/reference/excodeFormula.md);
  default `NULL`.

- weights_threshold:

  Numeric SD cutoff for
  [`surveillance::algo.farrington.assign.weights()`](https://rdrr.io/pkg/surveillance/man/algo.farrington.assign.weights.html)
  when down-weighting residuals in the baseline refit; default `2.58`.

- weights_threshold_baseline:

  Numeric Anscombe residual cutoff used to mark non-baseline points for
  clustering during initialization; default `1`.

- set_baseline_state:

  Logical; if `TRUE`, compute/assign the background state after
  initialization; default `FALSE`.

- maxIter:

  Integer maximum number of EM iterations; default `100`.

- verbose:

  Logical; print progress; default `FALSE`.

- return_full_model:

  Logical; if `FALSE` (default) return results only for the requested
  `timepoints`; if `TRUE`, return the full time series used for fitting
  at each `timepoint`.

## Value

A fitted
[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel-class.md)
containing, for the requested time points (or the full series if
`return_full_model = TRUE`), posterior state probabilities, fitted
means, p-values, Anscombe residuals, convergence information, and
information criteria.

## Details

The input `surv_ts` must be a `data.frame` with columns: **date**
(aggregation date), **observed** (case counts), optional **offset**
(e.g., population at risk), and optional **state** (e.g., 0 =
background, 1 = excess, `NA` = unknown). Passing `sts` objects (from
surveillance), lists of `sts`, or matrix-like collections of `sts` is
not supported. If the model’s emission uses a `MultiState`
`excodeFormula` that already contains a `surv_ts` slot, `run_excode()`
uses that data and, if missing, sets `timepoints` accordingly.

## See also

[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel-class.md),
[`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily-class.md),
[`excodeFormula`](https://robert-koch-institut.github.io/excode/reference/excodeFormula-class.md)

## Examples

``` r
data(shadar_df)
res_har_pois <- run_excode(
  surv_ts = shadar_df,
  timepoints = 295,
  distribution = "Poisson",
  states = 2,
  periodic_model = "Harmonic",
  time_trend = "Linear",
  set_baseline_state = TRUE
)
```
