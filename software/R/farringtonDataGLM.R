# The functions in this script are copied from the surveillance package https://cran.r-project.org/web/packages/surveillance/index.html as they are not exported in surveillance
# All functions in this script are needed to generate the data for the Farrington GLM
# The code was originally contributed to "surveillance" on 2013-06-26 by M. Salmon

algo.farrington.data.glm <- function(dayToConsider, b, freq,
                                     epochAsDate, epochStr,
                                     vectorOfDates, w, noPeriods,
                                     observed, population,
                                     verbose, pastWeeksNotIncluded, k) {
  # Identify reference time points

  # Same date but with one year, two year, etc, lag
  # b+1 because we need to have the current week in the vector
  referenceTimePoints <- algo.farrington.referencetimepoints(dayToConsider,
    b = b,
    freq = freq,
    epochAsDate = epochAsDate,
    epochStr = epochStr
  )

  if (!all(referenceTimePoints %in% vectorOfDates)) {
    ## previously only checked min(referenceTimePoints)
    stop(
      "Some reference time points did not exist; ",
      "decrease 'b' or postpone 'range'."
    )
  }

  if (verbose) {
    cat("k=", k, "\n")
  }

  # Create the blocks for the noPeriods between windows (including windows)
  # If noPeriods=1 this is a way of identifying windows, actually.

  blocks <- blocks(
    referenceTimePoints, vectorOfDates, epochStr, dayToConsider,
    b, w, noPeriods, epochAsDate
  )

  # Here add option for not taking the X past weeks into account
  # to avoid adaptation of the model to emerging outbreaks
  blocksID <- blocks
  blocksID[(k - pastWeeksNotIncluded):k] <- NA

  # Extract values for the timepoints of interest only

  blockIndexes <- which(is.na(blocksID) == FALSE)


  # Time

  # if epochAsDate make sure wtime has a 1 increment
  if (epochAsDate) {
    wtime <- (as.numeric(vectorOfDates[blockIndexes]) -
      as.numeric(vectorOfDates[blockIndexes][1])) / as.numeric(diff(vectorOfDates))[1]
  } else {
    wtime <- as.numeric(vectorOfDates[blockIndexes])
  }

  # Factors
  seasgroups <- as.factor(blocks[blockIndexes])

  # Observed
  response <- observed[blockIndexes]

  # Population
  pop <- population[blockIndexes]

  if (verbose) {
    print(response)
  }

  dataGLM <- data.frame(
    response = response, wtime = wtime, population = pop,
    seasgroups = seasgroups, vectorOfDates = vectorOfDates[blockIndexes]
  )
  dataGLM <- dataGLM[is.na(dataGLM$response) == FALSE, ]
  return(dataGLM)
}


algo.farrington.referencetimepoints <- function(dayToConsider, b = control$b, freq = freq, epochAsDate, epochStr) {
  if (epochAsDate) {
    referenceTimePoints <- as.Date(seq(
      as.Date(dayToConsider,
        origin = "1970-01-01"
      ),
      length.out = (b + 1), by = "-1 year"
    ))
  } else {
    referenceTimePoints <- seq(dayToConsider, length.out = (b + 1), by = -freq)
  }

  if (epochStr == "week") {
    # get the date of the Mondays/Tuesdays/etc so that it compares to
    # the reference data
    # (Mondays for Mondays for instance)

    # Vectors of same days near the date (usually the same week)
    # dayToGet
    dayToGet <- as.numeric(format(dayToConsider, "%w"))
    actualDay <- as.numeric(format(referenceTimePoints, "%w"))
    referenceTimePointsA <- referenceTimePoints -
      (actualDay
      - dayToGet)

    # Find the other "same day", which is either before or after referenceTimePoints

    referenceTimePointsB <- referenceTimePointsA + ifelse(referenceTimePointsA > referenceTimePoints, -7, 7)


    # For each year choose the closest Monday/Tuesday/etc
    # The order of referenceTimePoints is NOT important

    AB <- cbind(referenceTimePointsA, referenceTimePointsB)
    ABnumeric <- cbind(as.numeric(referenceTimePointsA), as.numeric(referenceTimePointsB))
    distMatrix <- abs(ABnumeric - as.numeric(referenceTimePoints))
    idx <- (distMatrix[, 1] > distMatrix[, 2]) + 1
    referenceTimePoints <- as.Date(AB[cbind(1:dim(AB)[1], idx)], origin = "1970-01-01")
  }

  return(referenceTimePoints)
}


blocks <- function(referenceTimePoints, vectorOfDates, freq, dayToConsider, b, w, p,
                   epochAsDate) {
  ## INPUT
  # freq: are we dealing with daily/weekly/monthly data?

  # b: how many years to go back in time

  # w: half window length around the reference timepoints

  # p: number of noPeriods one wants the year to be split into

  ## VECTOR OF ABSOLUTE NUMBERS
  # Very useful to write the code!

  vectorOfAbsoluteNumbers <- seq_len(length(vectorOfDates))


  # logical vector indicating where the referenceTimePoints
  # are in the vectorOfDates
  referenceTimePointsOrNot <- vectorOfDates %in% referenceTimePoints

  ## VECTOR OF FACTORS
  vectorOfFactors <- rep(NA_real_, length(vectorOfDates))

  ## SETTING THE FACTORS
  # Current week
  if (epochAsDate == FALSE) {
    now <- which(vectorOfDates == dayToConsider)
  } else {
    now <- which(vectorOfDates == as.Date(dayToConsider))
  }



  vectorOfFactors[(now - w):now] <- p

  # Reference weeks

  referenceWeeks <- rev(as.numeric(
    vectorOfAbsoluteNumbers[referenceTimePointsOrNot]
  ))

  for (i in 1:b) {
    # reference week

    refWeek <- referenceWeeks[i + 1]

    vectorOfFactors[(refWeek - w):(refWeek + w)] <- p

    # The rest is only useful if ones want factors, otherwise only have
    # reference timepoints like in the old algo.farrington

    if (p != 1) {
      # Number of time points to be shared between vectors
      period <- referenceWeeks[i] - 2 * w - 1 - refWeek

      # Check that p is not too big
      if (period < (p - (2 * w + 1))) {
        stop("Number of factors too big!")
      }

      # Look for the length of blocks

      lengthOfBlocks <- period %/% (p - 1)
      rest <- period %% (p - 1)

      vectorLengthOfBlocks <- rep(lengthOfBlocks, p - 1)

      # share the rest of the Euclidean division among the first blocks

      add <- seq_len(rest)
      vectorLengthOfBlocks[add] <- vectorLengthOfBlocks[add] + 1

      # slight transformation necessary for the upcoming code with cumsum
      vectorLengthOfBlocks <- c(0, vectorLengthOfBlocks)

      # fill the vector

      for (j in 1:(p - 1)) {
        vectorOfFactors[(refWeek + w + 1 + cumsum(vectorLengthOfBlocks)[j]):
        (refWeek + w + 1 + cumsum(vectorLengthOfBlocks)[j + 1] - 1)] <- j
      }
    }
  }

  ## DONE!

  return(vectorOfFactors) # indent
}
