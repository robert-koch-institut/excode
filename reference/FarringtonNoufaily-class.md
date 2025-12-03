# This class is a container for the parameterization of the FarringtonNoufaily models.

This class is a container for the parameterization of the
FarringtonNoufaily models.

## Slots

- `noPeriods`:

  Number of levels in the factor which creates bins in each year to
  model seasonal patterns.

- `w`:

  The number of weeks before and after the current week to include in
  the bin which contains the respective week in each year.

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
