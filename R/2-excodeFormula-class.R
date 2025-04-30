#' This class is a generic container for a formula which specifies the model used for excess count detection.
#'
#' @slot name Name of the model. Must be one of c("FarringtonNoufaily", "Harmonic", "Custom", "MultiState").
#' @slot shared_params Indicates whether model parameters are shared across multiple time series.
#' @slot coefficients Coefficients of a fitted model.
#' @slot failed_optim TRUE if an error occured during model fitting.
#'
#' @exportClass excodeFormula
setClass("excodeFormula",
  slots = c(
    name = "character",
    shared_params = "logical",
    coefficients = "numeric",
    failed_optim = "numeric"
  )
)


#' This class is a container for the parameterization of the FarringtonNoufaily models.
#'
#' @slot noPeriods Number of levels in the factor which creates bins in each year to model seasonal patterns.
#' @slot w The number of weeks before and after the current week to include in the bin which contains the respective week in each year.
#' @slot timeTrend Indicates whether a time trend should be included in the model.
#' @slot timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#' @slot offset TRUE if an offset should be included in the model.
#' @slot formula_bckg A formula which models the 'normal' (background) states.
#' @slot formula A formula which models which includes variable(s) to model 'excess' state(s).
#'
#' @exportClass excodeFormula
setClass("FarringtonNoufaily",
  contains = "excodeFormula",
  slots = c(
    noPeriods = "numeric",
    w = "numeric",
    timeTrend = "logical",
    timepoints_per_unit = "numeric",
    offset = "logical",
    formula_bckg = "character",
    formula = "character"
  ),
  prototype = list(
    name = "FarringtonNoufaily"
  )
)

#' Prints description of FarringtonNoufaily object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("FarringtonNoufaily"), function(object) {
  cat("excodeFormula: ", is(object)[1], "\n", sep = "")
})



#' This class is a container for the parameterization of the Harmonic models.
#'
#' @slot S Number of oscillations during one year.
#' @slot timeTrend Indicates whether a time trend should be included in the model.
#' @slot timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#' @slot offset TRUE if an offset should be included in the model.
#' @slot formula_bckg A formula which models the 'normal' (background) states.
#' @slot formula A formula which models which includes variable(s) to model 'excess' state(s).
#'
#' @exportClass excodeFormula
setClass("Harmonic",
  contains = "excodeFormula",
  slots = c(
    S = "numeric",
    timeTrend = "logical",
    timepoints_per_unit = "numeric",
    offset = "logical",
    formula_bckg = "character",
    formula = "character"
  ),
  prototype = list(
    name = "Harmonic"
  )
)


#' This class is a container for the parameterization of Spline models. This class is experimental and not recommended for use!
#'
#' @slot df_season Degrees of freedom for spline modeling seasonality. This class uses cubic splines, which are not recommended to model periodic patterns.
#' @slot df_trend Degrees of freedom for spline modeling the time trend
#' @slot timeTrend Indicates whether a time trend should be included in the model.
#' @slot timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#' @slot offset TRUE if an offset should be included in the model.
#' @slot formula_bckg A formula which models the 'normal' (background) states.
#' @slot formula A formula which models which includes variable(s) to model 'excess' state(s).
#'
#' @exportClass excodeFormula
setClass("Splines",
  contains = "excodeFormula",
  slots = c(
    df_season = "numeric",
    df_trend = "numeric",
    timepoints_per_unit = "numeric",
    offset = "logical",
    formula_bckg = "character",
    formula = "character"
  ),
  prototype = list(
    name = "Splines"
  )
)

#' Prints description of Harmonic object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("Harmonic"), function(object) {
  cat("excodeFormula: ", is(object)[1], "\n", sep = "")
})


#' Prints description of Splines object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("Splines"), function(object) {
  cat("excodeFormula: ", is(object)[1], "\n", sep = "")
})

#' This class is a container for the parameterization of the Mean models.
#'
#' @slot timeTrend Indicates whether a time trend should be included in the model.
#' @slot timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#' @slot offset TRUE if an offset should be included in the model.
#' @slot formula_bckg A formula which models the 'normal' (background) states.
#' @slot formula A formula which models which includes variable(s) to model 'excess' state(s).
#'
#' @exportClass excodeFormula
setClass("Mean",
  contains = "excodeFormula",
  slots = c(
    timeTrend = "logical",
    timepoints_per_unit = "numeric",
    offset = "logical",
    formula_bckg = "character",
    formula = "character"
  ),
  prototype = list(
    name = "Mean"
  )
)

#' Prints description of Mean object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("Mean"), function(object) {
  cat("excodeFormula: ", is(object)[1], "\n", sep = "")
})



#' This class is a container for the parameterization using external data.
#'
#' @slot data A data.frame containing variables which are used to model case counts. All variables in the data.frame will be used in the model. The data.frame has to have the same number of rows as time points in the time series.
#' @slot timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#' @slot offset TRUE if an offset should be included in the model.
#' @slot formula_bckg A formula which models the 'normal' (background) states.
#' @slot formula A formula which models which includes variable(s) to model 'excess' state(s).
#'
#'
#' @exportClass excodeFormula
setClass("Custom",
  contains = "excodeFormula",
  slots = c(
    data = "data.frame",
    timepoints_per_unit = "numeric",
    offset = "logical",
    formula_bckg = "character",
    formula = "character"
  ),
  prototype = list(
    name = "Custom"
  )
)


#' Prints description of Custom object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("Custom"), function(object) {
  cat("excodeFormula: ", is(object)[1], "\n", sep = "")
})




#' This class is a container for the parameterization of a MultiState model.
#'
#' @slot nStates The number of states of the model.
#' @slot data A data.frame containing variables which are used to model case counts. All variables in the data.frame will be used in the model. The data.frame has to have the same number of rows as time points in the time series.
#' @slot intercept TRUE if the model should include an intercept, FALSE otherwise.
#' @slot timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#' @slot offset TRUE if an offset should be included in the model.
#' @slot formula_bckg A formula which models the 'normal' (background) states.
#' @slot formula A formula which models which includes variable(s) to model 'excess' state(s).
#'
#'
#' @exportClass excodeFormula
setClass("MultiState",
  contains = "excodeFormula",
  slots = c(
    nStates = "numeric",
    data = "data.frame",
    intercept = "logical",
    timepoints_per_unit = "numeric",
    offset = "logical",
    formula_bckg = "character",
    formula = "character"
  ),
  prototype = list(
    name = "MultiState"
  )
)


#' Prints description of MultiState object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("MultiState"), function(object) {
  cat("excodeFormula: ", is(object)[1], "\n", sep = "")
})


#' @title Create a formula of the model for excess count detection
#'
#' @description
#' Constructs a model formula object for detecting excess counts using various model types,
#' such as Mean, Harmonic, Farrington-Noufaily, Custom, MultiState, and Splines.
#'
#' @param name Character. Name of the model. Must be one of c("Mean", "FarringtonNoufaily", "Harmonic", "Custom", "MultiState").
#' @param S Integer. Number of oscillations during one year for Harmonic models. Default is 1.
#' @param timepoints_per_unit Integer. Number of time points within the considered time unit (e.g. 52 for weekly observations in a year). Default is 52.
#' @param noPeriods Integer. Number of levels in the factor which creates bins in each year to model seasonal patterns. Only used in 'FarringtonNoufaily' models. Default is 10.
#' @param w Integer. The number of weeks before and after the current week to include in the bin which contains the respective week in each year. Only used in 'FarringtonNoufaily' models. Default is 3.
#' @param timeTrend Logical. Indicates whether a time trend should be included in the model. Used in 'Mean', 'Harmonic' and 'FarringtonNoufaily' models. Default is TRUE.
#' @param intercept Logical. TRUE if the model should include an intercept, FALSE otherwise. Only used with 'MultiState' models.Default is TRUE.
#' @param data data.frame. A data.frame containing variables which are used to model case counts. All variables in the data.frame will be used in the model. The data.frame has to have the same number of rows as time points in the time series. Default is NULL.
#' @param df_season Integer. Degrees of freedom for spline modeling seasonality. This class uses cubic splines, which are not recommended to model periodic patterns.  Only used for 'Splines' models. Default is 4.
#' @param df_trend Integer. Degrees of freedom for spline modeling the time trend. Only used for 'Splines' models. Default is 1.
#' @param nStates The number of states of the model. Default is NULL.
#' @param shared_params Logical. Indicates whether model parameters are shared across multiple time series. Default is FALSE.
#' @param offset Logical. Whether an offset should be included in the model. Default is FALSE.
#'
#' @return An S4 object of class \code{excodeFormula} representing the specified model formula.
#'
#' @seealso \code{\linkS4class{excodeFormula}}
#'
#' @examples
#'
#' excode_formula_har <- excodeFormula("Harmonic")
#'
#' @export
excodeFormula <- function(name,
                          S = 1, timepoints_per_unit = 52,
                          noPeriods = 10, w = 3,
                          timeTrend = TRUE,
                          intercept = TRUE,
                          data = NULL,
                          df_season = 4,
                          df_trend = 1,
                          nStates = NULL,
                          shared_params = FALSE,
                          offset = FALSE) {
  name <- as.character(name)
  failed_optim <- as.numeric(0)
  S <- as.numeric(S)
  timepoints_per_unit <- as.numeric(timepoints_per_unit)
  noPeriods <- as.numeric(noPeriods)
  w <- as.numeric(w)
  shared_params <- as.logical(shared_params)
  timeTrend <- as.logical(timeTrend)
  if (!is.null(data)) {
    data <- as.data.frame(data)
  } else if (name == "Custom" & is.null(data)) {
    stop("Must provide data.frame with excodeFormula: =", name, "\n", sep = "")
  } else if (name == "MultiState" & is.null(data)) {
    data <- data.frame()
  }

  obj <- NULL

  if (name == "Mean") {
    obj <- new(name,
      name = name,
      timeTrend = timeTrend, timepoints_per_unit = timepoints_per_unit,
      offset = offset, shared_params = shared_params,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "Harmonic") {
    obj <- new(name,
      name = name, S = S,
      timeTrend = timeTrend, timepoints_per_unit = timepoints_per_unit,
      offset = offset, shared_params = shared_params,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "FarringtonNoufaily") {
    obj <- new(name,
      name = name, noPeriods = noPeriods,
      w = w, timeTrend = timeTrend, timepoints_per_unit = timepoints_per_unit,
      offset = offset, shared_params = shared_params,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "Custom") {
    obj <- new(name,
      name = name, data = data,
      timepoints_per_unit = timepoints_per_unit, offset = offset, shared_params = shared_params,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "Splines") {
    obj <- new(name,
      name = name, df_season = df_season,
      df_trend = df_trend, timepoints_per_unit = timepoints_per_unit,
      offset = offset, shared_params = shared_params,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "MultiState") {
    obj <- new(name,
      name = name,
      nStates = as.numeric(nStates),
      intercept = intercept,
      data = data,
      timepoints_per_unit = timepoints_per_unit,
      offset = offset, shared_params = shared_params,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else {
    stop("name must be on of c(Mean, Harmonic, FarringtonNoufaily, Custom)\n")
  }

  obj
}
