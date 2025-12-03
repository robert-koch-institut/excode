# This class is a container for the parameterization of the Mean models.

This class is a container for the parameterization of the Mean models.

## Slots

- `timeTrend`:

  Indicates whether a time trend should be included in the model.

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
