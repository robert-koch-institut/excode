setGeneric("extractModelData", function(survts, model_struct,
                                        time_point_to_consider,
                                        time_units_back) {
  standardGeneric("extractModelData")
})



setMethod("extractModelData",
  signature = c(
    "sts",
    "FarringtonNoufaily",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit
    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * timepoints_per_unit + 1 + model_struct@w, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    observed <- survts@observed[, 1]
    state <- survts@state[, 1]
    epochAsDate <- survts@epochAsDate
    if (epochAsDate) {
      vectorOfDates <- as.Date(survts@epoch, origin = "1970-01-01")
    } else {
      vectorOfDates <- seq_len(length(observed))
    }
    dayToConsider <- vectorOfDates[time_point_to_consider]
    timepoints_per_unit <- survts@freq
    if (epochAsDate) {
      epochStr <- switch(as.character(timepoints_per_unit),
        "12" = "month",
        "52" = "week",
        "365" = "day"
      )
    } else {
      epochStr <- "none"
      allTimePoints <- allTimePoints[-1]
    }
    population <- surveillance::population(survts)

    # Create data for Farrington GLM
    modelData <- algo.farrington.data.glm(
      dayToConsider, time_units_back, timepoints_per_unit, epochAsDate,
      epochStr, vectorOfDates, model_struct@w, model_struct@noPeriods,
      observed, population, FALSE,
      pastWeeksNotIncluded = 0,
      time_point_to_consider
    )[, 1:4]

    modelData$rtime <- allTimePoints - 1
    modelData$true_state <- survts@state[modelData$rtime, 1]
    if (!epochAsDate) {
      modelData$wtime <- 0:(nrow(modelData) - 1)
    }

    currTimePointData <- data.frame(
      response = survts@observed[time_point_to_consider, 1],
      wtime = modelData$wtime[nrow(modelData)] + 1,
      population = survts@populationFrac[time_point_to_consider, 1],
      true_state = survts@state[time_point_to_consider, 1],
      seasgroups = model_struct@noPeriods,
      rtime = time_point_to_consider
    )
    modelData <- rbind(modelData, currTimePointData)
    modelData$date <- epoch(survts)[modelData$rtime]

    modelData
  }
)



setMethod("extractModelData",
  signature = c(
    "sts",
    "Mean",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit
    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * model_struct@timepoints_per_unit + 1, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    states <- survts@state[, 1]
    observed <- observed(survts)
    if (model_struct@offset) {
      offset_denom <- surveillance::population(survts)
    } else {
      offset_denom <- rep(1, length(observed))
    }
    modelData <- data.frame(
      response = observed[allTimePoints],
      wtime = 0:(length(allTimePoints) - 1),
      true_state = states[allTimePoints],
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )
    modelData$date <- epoch(survts)[allTimePoints]

    modelData
  }
)

setMethod("extractModelData",
  signature = c(
    "sts",
    "Harmonic",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit
    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * model_struct@timepoints_per_unit + 1, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    states <- survts@state[, 1]
    observed <- observed(survts)
    if (model_struct@offset) {
      offset_denom <- surveillance::population(survts)
    } else {
      offset_denom <- rep(1, length(observed))
    }
    modelData <- data.frame(
      response = observed[allTimePoints],
      wtime = 0:(length(allTimePoints) - 1),
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
            modelData[, sin_name] <- modelData[, sin_name] + sin(2 * pi * j * modelData$wtime / model_struct@timepoints_per_unit)
            modelData[, cos_name] <- modelData[, cos_name] + cos(2 * pi * j * modelData$wtime / model_struct@timepoints_per_unit)
          }
        }
      }
    }
    modelData$date <- epoch(survts)[allTimePoints]

    modelData
  }
)


setMethod("extractModelData",
  signature = c(
    "sts",
    "Splines",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    timepoints_per_unit <- model_struct@timepoints_per_unit
    allTimePoints <- rev(seq(time_point_to_consider, length = time_units_back * model_struct@timepoints_per_unit + 1, by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    states <- survts@state[, 1]
    observed <- observed(survts)
    if (model_struct@offset) {
      offset_denom <- surveillance::population(survts)
    } else {
      offset_denom <- rep(1, length(observed))
    }

    modelData <- data.frame(
      response = observed[allTimePoints],
      wtime = 0:(length(allTimePoints) - 1),
      true_state = states[allTimePoints],
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )



    if (model_struct@df_trend > 0) {
      trend_df <- data.frame(ns(
        1:nrow(modelData),
        model_struct@df_trend
      ))
      names(trend_df) <- paste0("wtime_", 1:model_struct@df_trend)
      modelData <- cbind(modelData, trend_df)
    }
    if (model_struct@df_season > 0) {
      season_df <- data.frame(ns(
        1:model_struct@timepoints_per_unit,
        model_struct@df_season
      ))
      len <- nrow(modelData)
      take <- rep(1:nrow(season_df), ceiling(len / model_struct@timepoints_per_unit))[1:nrow(modelData)]
      season_df <- season_df[take, ]
      names(season_df) <- paste0("season_", 1:model_struct@df_season)


      # season_df <- data.frame(cbind(sin(2*pi*modelData$wtime/model_struct@timepoints_per_unit),
      #                   cos(2*pi*modelData$wtime/model_struct@timepoints_per_unit)))
      # names(season_df) <- paste0("season_", 1:2)

      modelData <- cbind(modelData, season_df)
    }
    modelData$date <- epoch(survts)[allTimePoints]

    modelData
  }
)




setMethod("extractModelData",
  signature = c(
    "sts",
    "Custom",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    allTimePoints <- rev(seq(time_point_to_consider, length = round(time_units_back * model_struct@timepoints_per_unit + 1), by = -1))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    states <- survts@state[, 1]
    observed <- observed(survts)

    if (model_struct@offset) {
      offset_denom <- surveillance::population(survts)
    } else {
      offset_denom <- rep(1, length(observed))
    }

    modelData <- data.frame(
      response = observed[allTimePoints],
      true_state = states[allTimePoints],
      wtime = 0:(length(allTimePoints) - 1),
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )
    modelData$date <- epoch(survts)[allTimePoints]
    data <- model_struct@data[allTimePoints, , drop = F]
    modelData <- cbind(modelData, data)

    modelData
  }
)


setMethod("extractModelData",
  signature = c(
    "sts",
    "MultiState",
    "numeric", "numeric"
  ),
  function(survts, model_struct, time_point_to_consider, time_units_back) {
    allTimePoints <- rev(seq(time_point_to_consider,
      length = round(time_units_back * model_struct@timepoints_per_unit + 1),
      by = -1
    ))
    allTimePoints <- allTimePoints[allTimePoints > 0]

    states <- survts@state[, 1]
    observed <- observed(survts)

    if (model_struct@offset) {
      offset_denom <- surveillance::population(survts)
    } else {
      offset_denom <- rep(1, length(observed))
    }

    modelData <- data.frame(
      response = observed[allTimePoints],
      true_state = states[allTimePoints],
      rtime = allTimePoints,
      population = offset_denom[allTimePoints]
    )
    modelData$date <- epoch(survts)[allTimePoints]
    if (nrow(model_struct@data) > 0 & ncol(model_struct@data) > 0) {
      data <- model_struct@data[allTimePoints, , drop = F]
      modelData <- cbind(modelData, data)
    }

    rownames(modelData) <- NULL
    modelData$wtime <- 0:(nrow(modelData) - 1)
    modelData
  }
)




setGeneric("addDistrData", function(distribution, modelData) standardGeneric("addDistrData"))

setMethod("addDistrData",
  signature = c("Poisson", "data.frame"),
  function(distribution, modelData) {
    modelData
  }
)
setMethod("addDistrData",
  signature = c("NegBinom", "data.frame"),
  function(distribution, modelData) {
    if (distribution@shared_nb_size) {
      modelData$shared_nb_size <- "all"
    } else {
      modelData$shared_nb_size <- modelData$id
    }
    modelData
  }
)

setGeneric("prepareData", function(survts, hmm, time_point_to_consider,
                                   id, time_units_back = as.numeric(5),
                                   past_weeks_not_included_training = as.numeric(0),
                                   past_weeks_not_included_state = as.numeric(26),
                                   past_weeks_not_included_init = as.numeric(26)) {
  standardGeneric("prepareData")
})

setMethod("prepareData",
  signature = c(
    "sts", "excodeModel", "ANY",
    "ANY", "ANY", "ANY",
    "ANY", "ANY"
  ),
  function(survts, hmm, time_point_to_consider,
           id, time_units_back,
           past_weeks_not_included_training,
           past_weeks_not_included_state,
           past_weeks_not_included_init) {
    modelData <- extractModelData(
      survts, hmm@emission@excode_formula,
      time_point_to_consider, time_units_back
    )
    # modelData$denom = log(survts@populationFrac[modelData$rtime])
    modelData$id <- id
    modelData <- addDistrData(hmm@emission@distribution, modelData)

    modelData$init <- TRUE
    if (past_weeks_not_included_init > 0) {
      start_ind <- nrow(modelData) - past_weeks_not_included_init - 1
      modelData$init[start_ind:nrow(modelData)] <- FALSE
    }

    modelData$training <- TRUE
    if (past_weeks_not_included_training > 0) {
      start_ind <- nrow(modelData) - past_weeks_not_included_training
      modelData$training[start_ind:nrow(modelData)] <- FALSE
    }

    modelData$state_training <- TRUE
    if (past_weeks_not_included_state > 0) {
      start_ind <- nrow(modelData) - past_weeks_not_included_state - 1
      modelData$state_training[start_ind:nrow(modelData)] <- FALSE
    }

    modelData$curr_week <- FALSE
    modelData$curr_week[nrow(modelData)] <- TRUE

    modelData
  }
)

setMethod("prepareData",
  signature = c(
    "list", "excodeModel", "ANY",
    "ANY", "ANY", "ANY",
    "ANY", "ANY"
  ),
  function(survts, hmm, time_point_to_consider,
           id, time_units_back,
           past_weeks_not_included_training,
           past_weeks_not_included_state,
           past_weeks_not_included_init) {
    if (length(names(survts)) == 0) {
      names(survts) <- paste0("id", 1:length(survts))
    }

    modelData <- list()
    id_rep <- sapply(names(survts), function(x) rep(x, nrow(survts[[x]]@observed)))
    for (n in names(survts)) {
      curr_hmm <- hmm
      if ("data" %in% slotNames(hmm@emission@excode_formula)) {
        if (nrow(hmm@emission@excode_formula@data) > 0) {
          take_subset <- which(n == id_rep)
          curr_hmm@emission@excode_formula@data <- curr_hmm@emission@excode_formula@data[take_subset, , drop = FALSE]
        }
      }
      modelData[[n]] <- prepareData(survts[[n]], curr_hmm, time_point_to_consider,
        id = n, time_units_back,
        past_weeks_not_included_training,
        past_weeks_not_included_init
      )
      rownames(modelData[[n]]) <- NULL
    }

    modelData <- do.call("rbind", modelData)
    rownames(modelData) <- NULL
    modelData
  }
)
