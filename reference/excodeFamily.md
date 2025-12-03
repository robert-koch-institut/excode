# Create a family of probability distributions for excess count detection.

Create a family of probability distributions for excess count detection.

## Usage

``` r
excodeFamily(name, nb_size = NA)
```

## Arguments

- name:

  Name of the probability distribution that should be used. Either
  "Poisson" or "NegBinom".

- nb_size:

  Size parameter of the Negative Binomial distribution. Only relevant if
  'MultiState' model is used with a Negative Binomial dsitribution.

## Value

An
[`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily-class.md)
object.

## See also

[`excodeFamily`](https://robert-koch-institut.github.io/excode/reference/excodeFamily-class.md)
