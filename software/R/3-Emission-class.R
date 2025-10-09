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
    mu1 = "numeric",
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
      glm=empty_glm
    )
  } else {
    obj <- new(paste0("EmissionGLM", is(distribution)[1]),
      distribution = distribution,
      excode_formula = excode_formula,
      glm=empty_glm
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
            rA <- NULL
            mu <- model_data$mu[model_data$state==0]
            observed <- model_data$response[model_data$state==0]
            phi <- max(c(1, summary(emission@glm)$dispersion))
            numer <- 3 * (observed^(2/3) - mu^(2/3))
            denom <- 2 * sqrt(phi * mu^(1/3))
            rA <- numer / denom
            as.vector(rA)
          }
)


setMethod("zscores",
          signature = c("EmissionGLMNegBinom", "data.frame"),
          function(emission, model_data) {
            rA <- NULL
            
            glm_data <- model.frame(emission@glm)#model_data#
            glm_formula <- emission@excode_formula@formula
            col_pop <- which(names(glm_data)=="offset(log(population))")
            names(glm_data)[col_pop] <- "population"
            col_w <- which(names(glm_data)=="(weights)")
            names(glm_data)[col_w] <- "weights"
            
            glm_data$population <- exp(glm_data$population)
            glm_weights <- glm_data$weights
            glm_data <- glm_data[,-col_w]
            
            non_zero_w <- which(glm_weights>0)
            glm_qpois <- glm(as.formula(glm_formula), 
                             data=glm_data[non_zero_w,], 
                             weights=glm_weights[non_zero_w],
                             family=quasipoisson())
            
            mu <- predict(glm_qpois, newdata=model_data,
                          type="response")[model_data$state==0]
            observed <- model_data$response[model_data$state==0]
            phi <- max(c(1, summary(glm_qpois)$dispersion))
            numer <- 3 * (observed^(2/3) - mu^(2/3))
            denom <- 2 * sqrt(phi * mu^(1/3))
            rA <- numer / denom
            as.vector(rA)
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

