#' @title Emission Class
#'
#' @description The `Emission` class is a generic container for the model used in excess count detection. It defines the probability distributions that describe the emission process in the model.
#'
#' @slot name A character string giving the name of the emission.
#' @slot mu A numeric matrix where each column stores the mean (`mu`) for the respective hidden state of the model.
#' @slot excode_formula An object of class `excodeFormula` specifying the model formula used to describe the emission process.
#' @slot glm An object of class `glm` representing the fitted generalized linear model for the emission.
#'
#' @seealso [stats::glm()] for generalized linear models; [methods::setClass()] for S4 class definitions.
#'
#' @exportClass Emission
setClass("Emission",
  slots = c(
    name = "character",
    mu = "matrix",
    excode_formula = "excodeFormula",
    glm = "glm"
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
    distribution = "Poisson"
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
    distribution = "NegBinom"
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
  formula_bckg <- create_formula(distribution, excode_formula)
  formula <- paste0(formula_bckg, " + state")
  excode_formula@formula <- formula
  excode_formula@formula_bckg <- formula_bckg
  empty_glm <- list()
  class(empty_glm) <- "glm"

  obj <- NULL
  if (!is.null(initial_mu)) {
    obj <- new(paste0("EmissionGLM", is(distribution)[1]),
      distribution = distribution,
      excode_formula = excode_formula, mu = initial_mu,
      glm = empty_glm
    )
  } else {
    obj <- new(paste0("EmissionGLM", is(distribution)[1]),
      distribution = distribution,
      excode_formula = excode_formula,
      glm = empty_glm
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


#' @title Returns Anscombe residuals for a fitted Emission.
#'
#' @param emission An \code{\linkS4class{Emission}} object.
#' @returns A vector containing the Anscombe residuals for the fitted baseline state.
#'
#' @seealso \code{\linkS4class{Emission}}
#'
#' @keywords internal
#' @noRd
setGeneric("zscores", function(emission, model_data) standardGeneric("zscores"))

setMethod("zscores",
  signature = c("EmissionGLMPoisson", "data.frame"),
  function(emission, model_data) {
    glm_data <- model.frame(emission@glm) # model_data#
    col_w <- which(names(glm_data) == "(weights)")
    names(glm_data)[col_w] <- "weights"
    glm_weights <- glm_data$weights
    baseline_state <- which(model_data$state == 0)
    glm_weights <- glm_weights[baseline_state]
    formula_baseline <- as.formula(emission@excode_formula@formula_bckg)
    # print(c(length(gamma), sum(gamma)))
    downweighted <- (glm_weights < 1)

    gamma <- length(glm_weights) / (sum(glm_weights^(downweighted)))
    omega <- numeric(length(glm_weights))
    omega[downweighted] <- gamma * (glm_weights[downweighted])
    omega[!downweighted] <- gamma

    suppressWarnings(qpois_model <- glm(formula_baseline,
      data = model_data[baseline_state, ],
      family = quasipoisson(link = "log"),
      weights = omega
    ))
    phi <- max(summary(qpois_model)$dispersion, 1)
    qpois_anscombe_res <- surveillance::anscombe.residuals(qpois_model, phi)
    as.vector(qpois_anscombe_res)
  }
)


setMethod("zscores",
  signature = c("EmissionGLMNegBinom", "data.frame"),
  function(emission, model_data) {
    glm_data <- model.frame(emission@glm) # model_data#
    col_w <- which(names(glm_data) == "(weights)")
    names(glm_data)[col_w] <- "weights"
    glm_weights <- glm_data$weights
    baseline_state <- which(model_data$state == 0)
    glm_weights <- glm_weights[baseline_state]
    formula_baseline <- as.formula(emission@excode_formula@formula_bckg)
    # print(c(length(gamma), sum(gamma)))
    downweighted <- (glm_weights < 1)

    gamma <- length(glm_weights) / (sum(glm_weights^(downweighted)))
    omega <- numeric(length(glm_weights))
    omega[downweighted] <- gamma * (glm_weights[downweighted])
    omega[!downweighted] <- gamma

    suppressWarnings(qpois_model <- glm(formula_baseline,
      data = model_data[baseline_state, ],
      family = quasipoisson(link = "log"),
      weights = omega
    ))
    phi <- max(summary(qpois_model)$dispersion, 1)
    qpois_anscombe_res <- surveillance::anscombe.residuals(qpois_model, phi)
    as.vector(qpois_anscombe_res)
  }
)


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
    emission_df <- data.frame(emission@mu)
    names(emission_df) <- paste0("mu", 0:(ncol(emission_df) - 1))
    emission_summary <- emission_df
    if (!is.null(family_df)) {
      emission_summary <- cbind(emission_summary, family_df)
    }
    emission_summary
  }
)
