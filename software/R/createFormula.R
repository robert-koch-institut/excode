setGeneric("createFormula", function(distribution, excode_formula) standardGeneric("createFormula"))

setMethod("createFormula",
  signature = c("excodeFamily", "FarringtonNoufaily"),
  function(distribution, excode_formula) {
    params <- c("1", ifelse(excode_formula@timeTrend, "wtime", ""), "seasgroups")

    # Add shared model specification and collapse
    if (!excode_formula@shared_params) {
      params[params == "1"] <- "id"
      params[params != "id"] <- paste(params[params != "id"], "*id", sep = "")
    }
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
    wtime <- NULL
    if (excode_formula@timeTrend) {
      wtime <- "wtime"
    }
    params <- c("1", wtime, season)

    # Add shared model specification and collapse
    if (!excode_formula@shared_params) {
      params[params == "1"] <- "id"
      params[params != "id"] <- paste(params[params != "id"], "*id", sep = "")
    }
    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)



setMethod("createFormula",
  signature = c("excodeFamily", "Splines"),
  function(distribution, excode_formula) {
    season <- NULL
    if (excode_formula@df_season > 0) {
      season <- paste0("season_", 1:excode_formula@df_season)
    }
    wtime <- NULL
    if (excode_formula@df_trend > 0) {
      wtime <- paste0("wtime_", 1:excode_formula@df_trend)
    }
    params <- c("1", wtime, season)

    # Add shared model specification and collapse
    if (!excode_formula@shared_params) {
      params[params == "1"] <- "id"
      params[params != "id"] <- paste(params[params != "id"], "*id", sep = "")
    }
    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)

setMethod("createFormula",
  signature = c("excodeFamily", "Mean"),
  function(distribution, excode_formula) {
    wtime <- NULL
    if (excode_formula@timeTrend) {
      wtime <- "wtime"
    }
    params <- c("1", wtime)

    # Add shared model specification and collapse
    if (!excode_formula@shared_params) {
      params[params == "1"] <- "id"
      params[params != "id"] <- paste(params[params != "id"], "*id", sep = "")
    }
    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")

    paste0("response ~ ", params)
  }
)

setMethod("createFormula",
  signature = c("excodeFamily", "Custom"),
  function(distribution, excode_formula) {
    params <- c("1", names(excode_formula@data))

    # Add shared model specification and collapse
    if (!excode_formula@shared_params) {
      params[params == "1"] <- "id"
      params[params != "id"] <- paste(params[params != "id"], "*id", sep = "")
    }
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
    # Add shared model specification and collapse
    # if(!excode_formula@shared_params) {
    #  params[params=="1"] = "id"
    #  params[params!="id"] = paste(params[params!="id"], "*id", sep="")
    # }
    params[length(params) + 1] <- "offset(log(population))"
    params <- paste(params, collapse = " + ")
    paste0("response ~ ", params)
  }
)
