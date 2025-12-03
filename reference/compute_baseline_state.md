# Update EXCODE model by setting baseline states in surv_ts

Fits a quasi-Poisson GLM using the \*\*background formula\*\* stored in
`excode_multistate@emission@excode_formula@formula_bckg` on the data
found in `excode_multistate@emission@excode_formula@surv_ts` (optionally
augmented with covariates from `@data` if present), computes Anscombe
residuals, and \*\*writes\*\* a `state` column back into `@surv_ts`:
entries are `0` when residuals are below `weights_threshold_baseline`
(or the observed count is `0`), and `NA` otherwise. Returns the updated
`excodeModel`.

## Usage

``` r
compute_baseline_state(excode_multistate, weights_threshold_baseline = 1)
```

## Arguments

- excode_multistate:

  An `excodeModel` returned by
  [`init_excode()`](https://robert-koch-institut.github.io/excode/reference/init_excode.md)
  with a `"MultiState"` emission.

- weights_threshold_baseline:

  Numeric residual cutoff; residuals strictly below this threshold are
  labeled baseline (`0`). Default `1`.

## Value

The \*\*updated\*\* `excodeModel` with
`excode_multistate@emission@excode_formula@surv_ts$state` set.
