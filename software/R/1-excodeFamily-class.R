#' This class is a generic container for a family of probability distributions
#'
#' @slot name Name of the probability distribution.
#' @slot params Character vector of name of parameters.
#'
#' @exportClass excodeFamily
setClass("excodeFamily",
  slots = c(
    name = "character",
    params = "character"
  )
)


#' This class is a generic container for a family of Poisson distributions
#'
#' @exportClass Poisson
setClass("Poisson",
  contains = "excodeFamily",
  prototype = list(
    name = "Poisson",
    params = c("mu")
  )
)

#' Prints description of Poisson object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("Poisson"), function(object) {
  cat("excodeFamily: ", is(object)[1], "\n", sep = "")
})


#' This class is a generic container for a family of Negative Binomial distributions
#'
#' @slot nb_size Size parameter of the Negative Binomial distribution. Only relevant if 'MultiState' model is used with a Negative Binomial dsitribution.
#'
#' @exportClass NegBinom
setClass("NegBinom",
  contains = "excodeFamily",
  slots = c(
    nb_size = "numeric"
  ),
  prototype = list(
    name = "NegBinom",
    params = c("mu", "nb_size")
  )
)

#' Prints description of NegBinom object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("NegBinom"), function(object) {
  cat("excodeFamily: ", is(object)[1], "\n", sep = "")
})


#' @title Create a family of probability distributions for excess count detection.
#'
#' @param name Name of the probability distribution that should be used. Either "Poisson" or "NegBinom".
#' @param nb_size Size parameter of the Negative Binomial distribution. Only relevant if 'MultiState' model is used with a Negative Binomial dsitribution.
#' @returns An \code{\linkS4class{excodeFamily}} object.
#'
#' @seealso \code{\linkS4class{excodeFamily}}
#'
#' @export
excodeFamily <- function(name,
                         nb_size = NA) {
  name <- as.character(name)

  obj <- NULL

  if (name == "Poisson") {
    obj <- new(name, name = name)
  }
  if (name == "NegBinom") {
    obj <- new(name,
      name = name,
      nb_size = as.numeric(nb_size)
    )
  }

  obj
}


#' @title Returns parameter names of an excodeFamily object.
#'
#' @param excodeFamily An excodeFamily object.
#' @returns Character with the parameter names of the distribution.
#'
#' @seealso \code{\linkS4class{excodeFamily}}
#'
#' @keywords internal
#' @noRd
setGeneric("getParamNames", function(distribution) standardGeneric("getParamNames"))

setMethod("getParamNames",
  signature = c("excodeFamily"),
  function(distribution) {
    distribution@params
  }
)


#' @title Returns estimated parameters of an excodeFamily.
#'
#' @param family An excodeFamily object.
#' @returns NULL in case the excodeFamily is a Poisson distribution. The size parameter in case the excodeFamily is a Negative Binomial distribution.
#'
#' @seealso \code{\linkS4class{NegBinom}}, \code{\linkS4class{Poisson}}
#'
#' @keywords internal
#' @noRd
setGeneric("summary_family", function(family) standardGeneric("summary_family"))

setMethod("summary_family",
  signature = c("Poisson"),
  function(family) {
    NULL
  }
)

setMethod("summary_family",
  signature = c("NegBinom"),
  function(family) {
    data.frame(size = family@nb_size)
  }
)


#' @title Subsets size parameter of Negative Binomial distribution after model fitting to include the current timepoint only.
#'
#' @param excode_family An excodeFamily object.
#' @returns An excode_family object that contains size parameter for the current timepoint only.
#'
#' @seealso \code{\linkS4class{NegBinom}}, \code{\linkS4class{Poisson}}
#'
#' @keywords internal
#' @noRd
setGeneric("subset_nb_size", function(excode_family, ...) standardGeneric("subset_nb_size"))

setMethod("subset_nb_size",
  signature = c("NegBinom"),
  function(excode_family, index) {
    excode_family@nb_size <- excode_family@nb_size[index]
    excode_family
  }
)

setMethod("subset_nb_size",
  signature = c("Poisson"),
  function(excode_family, index) {
    excode_family
  }
)


#' @title Merge two excodeFamily after model fitting.
#'
#' @param excode_family1 An excodeFamily object.
#' @param excode_family1 An excodeFamily object.
#' @returns An excode_family object that contains size parameters of excode_family1 and excode_family2.
#'
#' @seealso \code{\linkS4class{NegBinom}}, \code{\linkS4class{Poisson}}
#'
#' @keywords internal
#' @noRd
setGeneric("merge_excode_family_result", function(excode_family1, excode_family2) standardGeneric("merge_excode_family_result"))

setMethod("merge_excode_family_result",
  signature = c("NegBinom"),
  function(excode_family1, excode_family2) {
    excode_family1@nb_size <- c(
      excode_family1@nb_size,
      excode_family2@nb_size
    )
    excode_family1
  }
)

setMethod("merge_excode_family_result",
  signature = c("Poisson"),
  function(excode_family1, excode_family2) {
    excode_family1
  }
)


#' @title Set the size parameter of a Negative Binomial distribution.
#'
#' @param excode_family An excodeFamily object.
#' @returns An excodeFamily object with the updated size parameter.
#'
#' @seealso \code{\linkS4class{NegBinom}}, \code{\linkS4class{Poisson}}
#'
#' @keywords internal
#' @noRd
setGeneric("set_nb_size", function(excode_family, ...) standardGeneric("set_nb_size"))

setMethod("set_nb_size",
  signature = c("NegBinom"),
  function(excode_family, model_outpout, index) {
    excode_family@nb_size <- model_outpout$nb_size[index]
    excode_family
  }
)

setMethod("set_nb_size",
  signature = c("Poisson"),
  function(excode_family, model_outpout, index) {
    excode_family
  }
)
