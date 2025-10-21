setGeneric("createFormula", function(distribution, excode_formula) standardGeneric("createFormula"))

setMethod("createFormula",
  signature = c("excodeFamily", "FarringtonNoufaily"),
  function(distribution, excode_formula) {
    params <- c("1", ifelse(excode_formula@timeTrend, "timepoint", ""), "seasgroups")

    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)


setMethod("createFormula",
  signature = c("excodeFamily", "Harmonic"),
  function(distribution, excode_formula) {
    season <- NULL
    if (excode_formula@S > 0) {
      season <- paste0(c("sin", "cos"), excode_formula@S)
    }
    timepoint <- NULL
    if (excode_formula@timeTrend) {
      timepoint <- "timepoint"
    }
    params <- c("1", timepoint, season)

    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)


setMethod("createFormula",
  signature = c("excodeFamily", "Mean"),
  function(distribution, excode_formula) {
    timepoint <- NULL
    if (excode_formula@timeTrend) {
      timepoint <- "timepoint"
    }
    params <- c("1", timepoint)

    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)

setMethod("createFormula",
  signature = c("excodeFamily", "Custom"),
  function(distribution, excode_formula) {
    params <- c("1", names(excode_formula@data))

    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)


setMethod("createFormula",
  signature = c("excodeFamily", "MultiState"),
  function(distribution, excode_formula) {
    offset_par <- ifelse(excode_formula@intercept, "1", "0")
    params <- c(offset_par, names(excode_formula@data))
    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")
    paste0("response ~ ", params)
  }
)

create_formula <- function(distribution, excode_formula) {
  params <- paste(excode_formula@params, collapse = " + ")
  paste0("response ~ ", params)
}
