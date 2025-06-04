setGeneric("calculatePvalue", function(distribution, result) standardGeneric("calculatePvalue"))

setMethod("calculatePvalue",
  signature = c("NegBinom", "data.frame"),
  function(distribution, result) {
    pnbinom(
      q = result$observed - 1,
      size = result$nb_size,
      mu = result$mu0,
      lower.tail = FALSE
    )
  }
)

setMethod("calculatePvalue",
  signature = c("Poisson", "data.frame"),
  function(distribution, result) {
    ppois(
      q = result$observed - 1,
      lambda = result$mu0,
      lower.tail = FALSE
    )
  }
)


setMethod("calculatePvalue",
  signature = c("NegBinom", "excodeModel"),
  function(distribution, result) {
    pnbinom(
      q = result@observed - 1,
      size = result@emission@distribution@nb_size,
      mu = result@emission@mu0,
      lower.tail = FALSE
    )
  }
)

setMethod("calculatePvalue",
  signature = c("Poisson", "excodeModel"),
  function(distribution, result) {
    ppois(
      q = result@observed - 1,
      lambda = result@emission@mu0,
      lower.tail = FALSE
    )
  }
)
