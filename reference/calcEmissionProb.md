# Calculate emission probabilities.

Calculate emission probabilities.

## Usage

``` r
calcEmissionProb(distribution, modelData)

# S4 method for class 'Poisson,data.frame'
calcEmissionProb(distribution, modelData)

# S4 method for class 'NegBinom,data.frame'
calcEmissionProb(distribution, modelData)
```

## Arguments

- distribution:

  A 'Poisson or 'NegBinom' object.

- modelData:

  Input data.

## Value

Emission probabilites as a matrix, where each column contains
probabilities of one state.

## See also

[`NegBinom`](https://robert-koch-institut.github.io/excode/reference/NegBinom-class.md),
[`Poisson`](https://robert-koch-institut.github.io/excode/reference/Poisson-class.md)
