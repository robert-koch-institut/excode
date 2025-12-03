# Initialize a multi-state EXCODE model from a surveillance time series

Builds an initial multi-state EXCODE model by fitting a baseline GLM to
the background component, computing Anscombe residuals with
Farrington-style down-weighting, clustering non-baseline residuals into
\`states - 1\` groups to seed latent states, and estimating
state-specific mean trajectories for \`initial_mu\`.

## Usage

``` r
init_excode(
  surv_ts,
  timepoint,
  time_units_back,
  distribution,
  states,
  time_trend = "None",
  periodic_model,
  period_length = 52,
  intercept = TRUE,
  covariate_df = NULL,
  weights_threshold = 2.58,
  weights_threshold_baseline = 1,
  set_baseline_state = FALSE
)
```

## Arguments

- surv_ts:

  data.frame of surveillance counts (response and any covariates
  required by the EXCODE formula); if an \`offset\` column is present it
  is used when building the MultiState formula.

- timepoint:

  Integer index indicating the right end of the training window passed
  to \`extractModelData()\`.

- time_units_back:

  Positive integer giving how many past time units (relative to
  \`timepoint\`) to include for initialization.

- distribution:

  Character scalar passed to \`excodeFamily()\` defining the emission
  distribution, e.g., \`"Poisson"\` or \`"NegBinom"\` (NB size inferred
  from the baseline fit).

- states:

  Integer (\>= 2) number of latent states; state \`0\` is treated as
  background.

- time_trend:

  Character time-trend setting forwarded to \`excodeFormula(name =
  periodic_model, time_trend = ...)\`, one of
  \`c("Linear","Spline1","Spline2","Spline3","Spline4","None")\`;
  default \`"None"\`.

- periodic_model:

  Character model name for the background formula, one of
  \`c("Mean","FarringtonNoufaily","Harmonic","Custom")\`.

- period_length:

  Integer mapped to \`timepoints_per_unit\` in \`excodeFormula\` (e.g.,
  \`52\` weekly points per year, \`7\` daily with weekly cycle, \`365\`
  daily with annual cycle); default 52.

- intercept:

  Logical forwarded to \`excodeFormula(name = "MultiState", intercept =
  ...)\`; default \`TRUE\`.

- covariate_df:

  data.frame or \`NULL\`; optional covariates used when constructing
  \`excodeFormula()\`; defaults to \`NULL\`.

- weights_threshold:

  Numeric SD cutoff used by
  \`surveillance::algo.farrington.assign.weights()\` to down-weight
  large residuals in the baseline refit; default \`2.58\`.

- weights_threshold_baseline:

  Numeric Anscombe residual threshold for identifying non-baseline time
  points that are clustered to seed the \`states - 1\` non-baseline
  states; default \`1\`.

- set_baseline_state:

  Logical; if \`TRUE\`, computes/assigns the background state after
  initialization; default \`FALSE\`.

## Value

An \`excodeModel\` whose emission is a \`MultiState\` formula with
\`nStates = states\`, initialized transition matrix and \`initial_mu\`
(and NB size when applicable).

## Details

Internally constructs a background \`excodeFormula(name =
periodic_model, ...)\`, fits a quasi-Poisson GLM with Farrington-style
weights, clusters residuals above \`weights_threshold_baseline\` via
\`kmeans()\` to initialize non-background states, then fits a GLM with a
\`state\` factor to obtain \`initial_mu\`; if \`distribution =
"NegBinom"\`, the NB size is initialized from the baseline fit.

## See also

[`excodeFormula`](https://robert-koch-institut.github.io/excode/reference/excodeFormula.md),
[`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily.md),
[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel.md),
[`surveillance::anscombe.residuals`](https://rdrr.io/pkg/surveillance/man/anscombe.residuals.html),
[`surveillance::algo.farrington.assign.weights`](https://rdrr.io/pkg/surveillance/man/algo.farrington.assign.weights.html)
