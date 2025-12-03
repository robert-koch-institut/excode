# Create a formula for excess count detection

Constructs an \`excodeFormula\` object for detecting excess counts using
model types such as Mean, Harmonic, Farrington-Noufaily, Custom, and
MultiState.

## Usage

``` r
excodeFormula(
  name,
  S = 1,
  timepoints_per_unit = 52,
  noPeriods = 10,
  w = 3,
  timeTrend = TRUE,
  time_trend = "Linear",
  intercept = TRUE,
  data = NULL,
  surv_ts = data.frame(),
  nStates = NULL,
  offset = FALSE
)
```

## Arguments

- name:

  Character; model name, one of
  \`c("Mean","FarringtonNoufaily","Harmonic","Custom","MultiState")\`.

- S:

  Integer; number of harmonic oscillations per year for Harmonic models;
  default 1.

- timepoints_per_unit:

  Integer; number of time points per unit (e.g., 52 weekly points per
  year, 365 daily points per year); default 52.

- noPeriods:

  Integer; number of seasonal bins per year (Farrington-Noufaily only);
  default 10.

- w:

  Integer; half-window (weeks) around the current week used to form
  seasonal bins (Farrington-Noufaily only); default 3.

- timeTrend:

  Logical; include a time trend (Mean, Harmonic, Farrington-Noufaily);
  default TRUE.

- time_trend:

  Character; time-trend form, one of
  \`c("Linear","Spline1","Spline2","None")\`; default "Linear".

- intercept:

  Logical; include an intercept (used in MultiState and Custom models);
  default TRUE.

- data:

  data.frame or NULL; covariate data used to build the model (all
  columns are included); required for "Custom", optional for others;
  defaults to NULL.

- surv_ts:

  data.frame; time-series data (used by MultiState); defaults to an
  empty data.frame.

- nStates:

  Integer or NULL; number of states for MultiState models; default NULL.

- offset:

  Logical; whether to include an offset (population) term; default
  FALSE.

## Value

An S4 object of class \`excodeFormula\` representing the specified model
formula.

## See also

[`excodeFormula`](https://robert-koch-institut.github.io/excode/reference/excodeFormula-class.md)
