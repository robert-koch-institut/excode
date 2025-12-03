# This class is a container for the parameterization of a MultiState model.

This class is a container for the parameterization of a MultiState
model.

## Slots

- `nStates`:

  The number of states of the model.

- `data`:

  A data.frame containing variables which are used to model case counts.
  All variables in the data.frame will be used in the model. The
  data.frame has to have the same number of rows as time points in the
  time series.

- `timepoints_per_unit`:

  Number of time points within the considered time unit (e.g. 52 for
  weekly observations in a year).

- `offset`:

  TRUE if an offset should be included in the model.

- `formula_bckg`:

  A formula which models the 'normal' (background) states.

- `formula`:

  A formula which models which includes variable(s) to model 'excess'
  state(s).
