#' @title excodeFormula Class
#' 
#' @description A generic container for a formula specifying the model used for excess count detection.
#' 
#' @slot name Character; name of the model, one of `c("FarringtonNoufaily","Harmonic","Custom","MultiState")`.
#' @slot params Character vector of parameter names included in the model formula.
#' @slot intercept Logical; `TRUE` if the model includes an intercept, `FALSE` otherwise.
#' @slot coefficients Numeric vector of fitted model coefficients.
#' @slot failed_optim Numeric indicator (0/1) showing whether an error occurred during model fitting.
#' @slot time_trend Character; type of time trend used in the model (e.g., `"Linear"`, `"Spline1"`, `"Spline2"`, `"None"`).
#' 
#' @exportClass excodeFormula
setClass("excodeFormula",
  slots = c(
    name = "character",
    params = "character",
    intercept = "logical",
    coefficients = "numeric",
    failed_optim = "numeric",
    time_trend = "character"
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


#' Prints description of Harmonic object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("Harmonic"), function(object) {
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
    surv_ts = "data.frame",
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


#' @title Create a formula for excess count detection
#'
#' @description Constructs an `excodeFormula` object for detecting excess counts using model types such as Mean, Harmonic, Farrington-Noufaily, Custom, and MultiState.
#'
#' @param name Character; model name, one of `c("Mean","FarringtonNoufaily","Harmonic","Custom","MultiState")`.
#' @param S Integer; number of harmonic oscillations per year for Harmonic models; default 1.
#' @param timepoints_per_unit Integer; number of time points per unit (e.g., 52 weekly points per year, 365 daily points per year); default 52.
#' @param noPeriods Integer; number of seasonal bins per year (Farrington-Noufaily only); default 10.
#' @param w Integer; half-window (weeks) around the current week used to form seasonal bins (Farrington-Noufaily only); default 3.
#' @param timeTrend Logical; include a time trend (Mean, Harmonic, Farrington-Noufaily); default TRUE.
#' @param time_trend Character; time-trend form, one of `c("Linear","Spline1","Spline2","None")`; default "Linear".
#' @param intercept Logical; include an intercept (used in MultiState and Custom models); default TRUE.
#' @param data data.frame or NULL; covariate data used to build the model (all columns are included); required for "Custom", optional for others; defaults to NULL.
#' @param surv_ts data.frame; time-series data (used by MultiState); defaults to an empty data.frame.
#' @param nStates Integer or NULL; number of states for MultiState models; default NULL.
#' @param offset Logical; whether to include an offset (population) term; default FALSE.
#'
#' @return An S4 object of class `excodeFormula` representing the specified model formula.
#'
#' @seealso \code{\linkS4class{excodeFormula}}
#'
#' @export
excodeFormula <- function(name,
                          S = 1, timepoints_per_unit = 52,
                          noPeriods = 10, w = 3,
                          timeTrend = TRUE,
                          time_trend = "Linear",
                          intercept = TRUE,
                          data = NULL,
                          surv_ts = data.frame(),
                          nStates = NULL,
                          offset = FALSE) {
  name <- as.character(name)
  failed_optim <- as.numeric(0)
  S <- as.numeric(S)
  timepoints_per_unit <- as.numeric(timepoints_per_unit)
  noPeriods <- as.numeric(noPeriods)
  w <- as.numeric(w)
  timeTrend <- as.logical(timeTrend)
  time_trend <- as.character(time_trend)
  if (!is.null(data)) {
    data <- as.data.frame(data)
  } else if (name == "Custom" & is.null(data)) {
    stop("Must provide data.frame with excodeFormula: =", name, "\n", sep = "")
  } else if (name == "MultiState" & is.null(data)) {
    data <- data.frame()
  }

  obj <- NULL
  spline_params <- ""
  time_trend_orig <- time_trend
  if (length(grep("Spline", time_trend)) > 0) {
    n_knots <- as.numeric(gsub("Spline", "", time_trend))
    spline_params <- paste0("t_spline", 1:(n_knots + 1))
    time_trend <- "Spline"
  }
  params_time <- switch(time_trend,
    Linear = "timepoint",
    Spline = spline_params,
    None = NULL
  )
  time_trend <- time_trend_orig
  params <- c(ifelse(intercept, "1", "0"), params_time)
  if (name == "Mean") {
    params <- c(params, "offset(log(population))")
    obj <- new(name,
      name = name,
      timeTrend = timeTrend,
      time_trend = time_trend,
      params = params,
      intercept = intercept,
      timepoints_per_unit = timepoints_per_unit,
      offset = offset,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "Harmonic") {
    params <- c(
      params,
      paste0("sin", 1:S), paste0("cos", 1:S),
      "offset(log(population))"
    )
    obj <- new(name,
      name = name, S = S,
      timeTrend = timeTrend,
      time_trend = time_trend,
      params = params,
      intercept = intercept,
      timepoints_per_unit = timepoints_per_unit,
      offset = offset,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "FarringtonNoufaily") {
    params <- c(
      params, "seasgroups",
      "offset(log(population))"
    )
    obj <- new(name,
      name = name, noPeriods = noPeriods,
      w = w, timeTrend = timeTrend,
      time_trend = time_trend,
      params = params,
      intercept = intercept,
      timepoints_per_unit = timepoints_per_unit,
      offset = offset,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "Custom") {
    params <- c(colnames(data), "offset(log(population))")
    obj <- new(name,
      name = name, data = data,
      params = params,
      intercept = intercept,
      timepoints_per_unit = timepoints_per_unit, offset = offset,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else if (name == "MultiState") {
    params <- c(colnames(data), "offset(log(population))")
    obj <- new(name,
      name = name,
      nStates = as.numeric(nStates),
      intercept = intercept,
      data = data,
      surv_ts = surv_ts,
      params = params,
      timepoints_per_unit = timepoints_per_unit,
      offset = offset,
      coefficients = as.numeric(NA),
      failed_optim = failed_optim
    )
  } else {
    stop("name must be on of c(Mean, Harmonic, FarringtonNoufaily, Custom, MultiState)\n")
  }

  obj
}
