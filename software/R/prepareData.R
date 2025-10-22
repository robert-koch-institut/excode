setGeneric("extractModelData", function(survts, model_struct,
                                        time_point_to_consider,
                                        time_units_back) {
  standardGeneric("extractModelData")
})


setMethod("extractModelData",
  signature = c(
    "data.frame",
    "FarringtonNoufaily",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit

    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * timepoints_per_unit + 1 + model_struct@w, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    observed <- survts$observed
    if (!"state" %in% names(survts)) {
      survts$state <- NA
    }
    state <- survts$state
    vectorOfDates <- survts$date
    dayToConsider <- vectorOfDates[time_point_to_consider]

    epochStr <- switch(as.character(timepoints_per_unit),
      "12" = "month",
      "52" = "week",
      "365" = "day"
    )
    if (model_struct@offset) {
      offset_denom <- survts$offset
    } else {
      offset_denom <- rep(1, length(observed))
    }
    survts$offset <- offset_denom
    population <- survts$offset
    # Create data for Farrington GLM

    modelData <- algo.farrington.data.glm(
      dayToConsider, time_units_back, timepoints_per_unit, TRUE,
      epochStr, vectorOfDates, model_struct@w, model_struct@noPeriods,
      observed, population, FALSE,
      pastWeeksNotIncluded = 0,
      time_point_to_consider
    )[, 1:4]

    modelData$rtime <- allTimePoints - 1
    modelData$true_state <- survts$state[modelData$rtime]

    names(modelData)[names(modelData) == "wtime"] <- "timepoint"

    currTimePointData <- data.frame(
      response = survts$observed[time_point_to_consider],
      timepoint = modelData$timepoint[nrow(modelData)] + 1,
      population = survts$offset[time_point_to_consider],
      true_state = survts$state[time_point_to_consider],
      seasgroups = model_struct@noPeriods,
      rtime = time_point_to_consider
    )
    modelData <- rbind(modelData, currTimePointData)
    modelData$date <- survts$date[modelData$rtime]
    modelData <- add_splines(modelData, model_struct@time_trend)

    modelData
  }
)


setMethod("extractModelData",
  signature = c(
    "data.frame",
    "Mean",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit
    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * model_struct@timepoints_per_unit + 1, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]
    if (!"state" %in% names(survts)) {
      survts$state <- NA
    }
    states <- survts$state
    observed <- survts$observed
    if (model_struct@offset) {
      offset_denom <- survts$offset
    } else {
      offset_denom <- rep(1, length(observed))
    }
    modelData <- data.frame(
      response = observed[allTimePoints],
      timepoint = 0:(length(allTimePoints) - 1),
      true_state = states[allTimePoints],
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )
    modelData$date <- survts$date[allTimePoints]
    modelData <- add_splines(modelData, model_struct@time_trend)

    modelData
  }
)

setMethod("extractModelData",
  signature = c(
    "data.frame",
    "Harmonic",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit
    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * model_struct@timepoints_per_unit + 1, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    if (!"state" %in% names(survts)) {
      survts$state <- NA
    }
    states <- survts$state
    observed <- survts$observed
    if (model_struct@offset) {
      offset_denom <- survts$offset
    } else {
      offset_denom <- rep(1, length(observed))
    }
    modelData <- data.frame(
      response = observed[allTimePoints],
      timepoint = 0:(length(allTimePoints) - 1),
      true_state = states[allTimePoints],
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )

    if (model_struct@S > 0) {
      for (i in 1:model_struct@S) {
        sin_name <- paste("sin", i, sep = "")
        cos_name <- paste("cos", i, sep = "")
        modelData[, sin_name] <- 0
        modelData[, cos_name] <- 0
        for (j in 1:i) {
          if (j == model_struct@S) {
            modelData[, sin_name] <- modelData[, sin_name] + sin(2 * pi * j * modelData$timepoint / model_struct@timepoints_per_unit)
            modelData[, cos_name] <- modelData[, cos_name] + cos(2 * pi * j * modelData$timepoint / model_struct@timepoints_per_unit)
          }
        }
      }
    }
    modelData$date <- survts$date[allTimePoints]
    modelData <- add_splines(modelData, model_struct@time_trend)

    modelData
  }
)


setMethod("extractModelData",
  signature = c(
    "data.frame",
    "Custom",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    allTimePoints <- rev(seq(time_point_to_consider, length = round(time_units_back * model_struct@timepoints_per_unit + 1), by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]
    if (!"state" %in% names(survts)) {
      survts$state <- NA
    }
    states <- survts$state
    observed <- survts$observed

    if (model_struct@offset) {
      offset_denom <- survts$offset
    } else {
      offset_denom <- rep(1, length(observed))
    }

    modelData <- data.frame(
      response = observed[allTimePoints],
      true_state = states[allTimePoints],
      timepoint = 0:(length(allTimePoints) - 1),
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )
    modelData$date <- survts$date[allTimePoints]
    data <- model_struct@data[allTimePoints, , drop = F]
    modelData <- cbind(modelData, data)

    modelData
  }
)


setMethod("extractModelData",
  signature = c(
    "data.frame",
    "MultiState",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    allTimePoints <- rev(seq(time_point_to_consider,
      length = round(time_units_back * model_struct@timepoints_per_unit + 1),
      by = -1
    ))
    allTimePoints <- 1:nrow(survts)
    allTimePoints <- allTimePoints[allTimePoints > 0]
    if (!"state" %in% names(survts)) {
      survts$state <- NA
    }
    states <- survts$state
    observed <- survts$observed

    if (model_struct@offset) {
      offset_denom <- survts$offset
    } else {
      offset_denom <- rep(1, length(observed))
    }

    modelData <- data.frame(
      response = observed[allTimePoints],
      true_state = states[allTimePoints],
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )
    modelData$date <- survts$date[allTimePoints]
    if (nrow(model_struct@data) > 0 & ncol(model_struct@data) > 0) {
      data <- model_struct@data[allTimePoints, , drop = F]
      modelData <- cbind(modelData, data)
    }

    rownames(modelData) <- NULL
    modelData$timepoint <- 0:(nrow(modelData) - 1)

    modelData
  }
)


setGeneric("prepareData", function(survts, hmm, time_point_to_consider,
                                   time_units_back = as.numeric(5)) {
  standardGeneric("prepareData")
})

setMethod("prepareData",
  signature = c(
    "data.frame", "excodeModel",
    "ANY", "ANY"
  ),
  function(survts, hmm, time_point_to_consider,
           time_units_back) {
    modelData <- extractModelData(
      survts, hmm@emission@excode_formula,
      time_point_to_consider, time_units_back
    )

    modelData$curr_week <- FALSE
    modelData$curr_week[nrow(modelData)] <- TRUE

    modelData
  }
)


add_splines <- function(model_data, time_trend) {
  # if (!time_trend %in% c("Linear", "Spline1", "Spline2", "None")) {
  #  stop("Invalid value for 'time_trend'. Must be one of: 'Linear', 'Spline1', 'Spline2', 'None'.")
  # }

  timepoint <- model_data$timepoint
  if (length(grep("Spline", time_trend)) == 1) { # %in% c("Spline1", "Spline2")) {
    n_knots <- as.numeric(gsub("Spline", "", time_trend))
    t_knots <- round(seq(1, length(timepoint),
      length = n_knots + 2
    )[-c(1, n_knots + 2)])
    spline_df <- as.data.frame(ns(timepoint,
      knots = t_knots
    ))
    names(spline_df) <- paste0("t_spline", 1:ncol(spline_df))
    spline_df$timepoint <- timepoint
    model_data <- dplyr::left_join(model_data,
      spline_df,
      by = "timepoint"
    )
  }
  model_data
}
