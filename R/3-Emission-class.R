#' This class is a generic container for the model used for excess count detection, which defines the probability distributions.
#'
#' @slot name Name of the emission.
#' @slot mu Matrix that stores in each column the mean for the respective state.
#' @slot mu0 Numeric vector containing the mean for the 'normal' state in 2-state models.
#' @slot mu1 Numeric vector containing the mean for the 'excess' state in 2-state models. Name of the emission.
#'
#' @exportClass Emission
setClass("Emission",
  slots = c(
    name = "character",
    mu = "matrix",
    mu0 = "numeric",
    mu1 = "numeric"
  )
)


#' This class defines the Emission function for Poisson GLMs.
#'
#' @slot distribution Name of the distribution ("Poisson").
#' @slot excode_formula An \code{\linkS4class{excodeFormula}} object.
#'
#' @seealso \code{\linkS4class{excodeFormula}}
#'
#' @exportClass EmissionGLMPoisson
setClass("EmissionGLMPoisson",
  contains = "Emission",
  slots = c(
    distribution = "Poisson",
    excode_formula = "excodeFormula"
  )
)


#' This class defines the Emission function for Negative Binomial GLMs.
#'
#' @slot distribution Name of the distribution ("NegBinom").
#' @slot excode_formula An \code{\linkS4class{excodeFormula}} object.
#'
#' @seealso \code{\linkS4class{excodeFormula}}
#'
#' @exportClass EmissionGLMNegBinom
setClass("EmissionGLMNegBinom",
  contains = "Emission",
  slots = c(
    distribution = "NegBinom",
    excode_formula = "excodeFormula"
  )
)


#' @title Create an Emission object for excess count detection
#'
#' @param distribution excodeFamily object.
#' @param excode_formula excodeFormula object.
#' @param initial_mu Initial estimates of the mean for 'MultiState' models.
#'
#' @returns An \code{\linkS4class{Emission}} object.
#'
#' @export
Emission <- function(distribution, excode_formula, initial_mu = NULL) {
  formula_bckg <- createFormula(distribution, excode_formula)
  formula <- paste0(formula_bckg, " + state")
  excode_formula@formula <- formula
  excode_formula@formula_bckg <- formula_bckg

  obj <- NULL
  if (!is.null(initial_mu)) {
    obj <- new(paste0("EmissionGLM", is(distribution)[1]),
      distribution = distribution,
      excode_formula = excode_formula, mu = initial_mu
    )
  } else {
    obj <- new(paste0("EmissionGLM", is(distribution)[1]),
      distribution = distribution,
      excode_formula = excode_formula
    )
  }
}


#' @title Update Emission function during model fitting
#'
#' @param emission An \code{\linkS4class{Emission}} object.
#' @param dat A data.frame containing the data for fitting.
#' @param omega Regression weights (posterior probabilities).
#'
#' @returns#' @returns A list containing \itemize{
#'    \item{model - }{The data.}
#'    \item{emission - }{The updated emission}
#' }
#'
#' @seealso \code{\linkS4class{Emission}}
#'
#' @keywords internal
#' @noRd
setGeneric("updateEmission", function(emission, dat, omega) standardGeneric("updateEmission"))



#' @title Returns estimated parameters of an excodeFamily.
#'
#' @param emission An \code{\linkS4class{Emission}} object.
#' @returns A data.frame containing the summary of the Emission object.
#'
#' @seealso \code{\linkS4class{Emission}}
#'
#' @keywords internal
#' @noRd
setGeneric("summary_emission", function(emission) standardGeneric("summary_emission"))

setMethod("summary_emission",
  signature = c("Emission"),
  function(emission) {
    family_df <- summary_family(emission@distribution)
    # emission_df <- data.frame(mu0=emission@mu0,
    #                          mu1=emission@mu1)
    emission_df <- data.frame(emission@mu)
    names(emission_df) <- paste0("mu", 0:(ncol(emission_df) - 1))
    emission_summary <- emission_df
    if (!is.null(family_df)) {
      emission_summary <- cbind(emission_summary, family_df)
    }
    emission_summary
  }
)
