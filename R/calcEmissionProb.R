#' @title Calculate emission probabilities.
#'
#' @param distribution A 'Poisson or 'NegBinom' object.
#' @param modelData Input data.
#'
#' @returns Emission probabilites as a matrix, where each column contains probabilities of one state.
#'
#' @seealso \code{\linkS4class{NegBinom}}, \code{\linkS4class{Poisson}}
#' @rdname calcEmissionProb
#'
#' @keywords internal
setGeneric("calcEmissionProb", function(distribution, modelData) standardGeneric("calcEmissionProb"))

#' @rdname calcEmissionProb
setMethod("calcEmissionProb",
  signature = c("Poisson", "data.frame"),
  function(distribution, modelData) {
    prob <- dpois(modelData$response, lambda = modelData$mu)
    emissionProb <- matrix(prob, ncol = length(unique(modelData$state)))
    emissionProb
  }
)

#' @rdname calcEmissionProb
setMethod("calcEmissionProb",
  signature = c("NegBinom", "data.frame"),
  function(distribution, modelData) {
    prob <- dnbinom(modelData$response, mu = modelData$mu, size = modelData$nb_size)
    emissionProb <- matrix(prob, ncol = length(unique(modelData$state)))
    emissionProb
  }
)
