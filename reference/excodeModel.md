# Create a model for excess count detection

Creates an object of class `excodeModel` for modeling excess counts
using a Hidden Markov Model (HMM) with user-defined emission
distributions and formula for modeling observed counts.

## Usage

``` r
excodeModel(
  family,
  formula,
  initial_mu = NULL,
  transMat = NULL,
  initProb = NULL,
  transMat_prior = TRUE,
  setBckgState = TRUE
)
```

## Arguments

- family:

  An
  [`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily-class.md)
  object defining the emission distribution (e.g., Poisson, Negative
  Binomial).

- formula:

  An
  [`excodeFormula`](https://robert-koch-institut.github.io/excode/reference/excodeFormula-class.md)
  object specifying the structure of the model (e.g., time trends,
  seasonality, ...).

- initial_mu:

  Initial estimates of the mean for 'MultiState' models.

- transMat:

  Inital transition probabilities.

- initProb:

  A numeric vector containing initial state probabilities (probabilities
  of of states at first time point) of the hidden Markov model.

- transMat_prior:

  Logical. Should a prior distribution be used for estimating transition
  probabilities? Default is `TRUE`.

- setBckgState:

  Logical. Should a background state be inferred for model fitting?
  Background states are initialized to 0 for time points with Anscombe
  residuals \< 1 from an initial model. Default is `TRUE`.

## Value

An object of class
[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel-class.md).

## See also

[`excodeModel`](https://robert-koch-institut.github.io/excode/reference/excodeModel-class.md),
[`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily-class.md),
[`excodeFormula`](https://robert-koch-institut.github.io/excode/reference/excodeFormula-class.md)
