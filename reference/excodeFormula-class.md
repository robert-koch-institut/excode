# excodeFormula Class

A generic container for a formula specifying the model used for excess
count detection.

## Slots

- `name`:

  Character; name of the model, one of
  \`c("FarringtonNoufaily","Harmonic","Custom","MultiState")\`.

- `params`:

  Character vector of parameter names included in the model formula.

- `intercept`:

  Logical; \`TRUE\` if the model includes an intercept, \`FALSE\`
  otherwise.

- `coefficients`:

  Numeric vector of fitted model coefficients.

- `failed_optim`:

  Numeric indicator (0/1) showing whether an error occurred during model
  fitting.

- `time_trend`:

  Character; type of time trend used in the model (e.g., \`"Linear"\`,
  \`"Spline1"\`, \`"Spline2"\`, \`"None"\`).
