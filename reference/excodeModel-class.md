# excodeModel Class

The \`excodeModel\` class is a generic container for the model used for
outbreak detection based on a hidden Markov model framework.

## Slots

- `nStates`:

  Number of states in the hidden Markov model.

- `emission`:

  Emission component of the hidden Markov model, including an
  \`excodeFamily\` and an \`excodeFormula\` object.

- `transitions`:

  Transition probabilities of the hidden Markov model.

- `transitions_prior`:

  Logical indicating whether a prior for transition probabilities should
  be used.

- `prior_weights`:

  Dirichlet prior weights for transition probabilities; cannot be
  changed by the user.

- `loglik_transitions`:

  Prior log-likelihood of transition probabilities; only relevant for
  model fitting.

- `initial_prob`:

  Initial state probabilities of the hidden Markov model.

- `setBckgState`:

  Logical indicating whether a background state should be computed for
  model fitting; background states are set to 0 for time points where
  the Anscombe residuals of the initial model are \< 1.

- `converged`:

  Logical indicating whether the EM algorithm converged.

- `niter`:

  Number of iterations used for model fitting.

- `LogLik`:

  Log-likelihood of the fitted model.

- `AIC`:

  Akaike information criterion of the fitted model.

- `BIC`:

  Bayesian information criterion of the fitted model.

- `posterior`:

  Matrix of posterior probabilities for each time point (rows) and state
  (columns).

- `pval`:

  P-values for testing whether each time point is in a normal or excess
  state.

- `zscore`:

  Standardized Anscombe residuals (z-scores) for each time point.

- `timepoint_fit`:

  Numeric time point in the time series used for fitting.

- `timepoint`:

  Numeric time point in the time series.

- `date`:

  Date corresponding to each time point.

- `observed`:

  Vector of observed counts or case numbers.

- `population`:

  Population size; default is 1.

- `error`:

  Character string containing an error message, if an error occurred
  during model fitting.
