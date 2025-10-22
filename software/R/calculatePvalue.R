setGeneric("calculatePvalue", function(distribution, result) standardGeneric("calculatePvalue"))


setMethod("calculatePvalue",
  signature = c("NegBinom", "excodeModel"),
  function(distribution, result) {
    pnbinom(
      q = result@observed - 1,
      size = result@emission@distribution@nb_size,
      mu = result@emission@mu[, 1],
      lower.tail = FALSE
    )
  }
)

setMethod("calculatePvalue",
  signature = c("Poisson", "excodeModel"),
  function(distribution, result) {
    pval <- ppois(
      q = result@observed - 1,
      lambda = result@emission@mu[, 1],
      lower.tail = FALSE
    )

    # size=mu/(phi-1)
    phi <- max(c(1, summary(result@emission@glm)$dispersion))
    if (phi > 1) {
      pval <- pnbinom(
        q = result@observed - 1,
        mu = result@emission@mu[, 1],
        size = result@emission@mu[, 1] / (phi - 1),
        lower.tail = FALSE
      )
    }
    pval
  }
)
