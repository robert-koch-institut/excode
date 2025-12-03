# Emission Class

The \`Emission\` class is a generic container for the model used in
excess count detection. It defines the probability distributions that
describe the emission process in the model.

## Slots

- `name`:

  A character string giving the name of the emission.

- `mu`:

  A numeric matrix where each column stores the mean (\`mu\`) for the
  respective hidden state of the model.

- `excode_formula`:

  An object of class \`excodeFormula\` specifying the model formula used
  to describe the emission process.

- `glm`:

  An object of class \`glm\` representing the fitted generalized linear
  model for the emission.

## See also

\[stats::glm()\] for generalized linear models; \[methods::setClass()\]
for S4 class definitions.
