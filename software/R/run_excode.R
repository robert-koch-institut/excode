#' @title Detect excess counts from single or multiple epidemiological time series. Internal function which is applied to a surveillance sts time series.
#'
#' @param surv_ts A surveillance time series (sts) object created using the surveillance package.
#' @param excode_model An object of class \code{\linkS4class{excodeModel}} which specifies the model parameters and structure
#' @param timepoints integer or sequence of integers which specifies the time points for which excess count detection should be performed
#' @param learning_type Indicates the type of learning, one of c("unsupervised", "semisupervised", "supervised")
#' @param maxIter Maximal number of iteration for EM-algorithm.
#' @param verbose Logical indicating wehther progress should be printed.
#' @param return_full_model If FALSE (default) only results for 'timepoint' are returend. If TRUE, the complete time series used for model fitting is returned by the function.
#' @param past_timepoints_not_incuded Past time points not included in initialization.
#' @param time_units_back Number of years to be used for model fitting.
#'
#' @seealso \code{\linkS4class{excodeModel}}
#'
#' @keywords internal
#' @noRd
run_excode_internal <- function(surv_ts, excode_model, timepoints, learning_type = "unsupervised", maxIter = 100,
                                verbose = FALSE, return_full_model = FALSE, past_timepoints_not_incuded = 26, time_units_back = 5) {
  excode_model_init <- excode_model

  transMat_init <- excode_model@transitions
  initProb_init <- excode_model@initial_prob

  all_errors <- rep("", length(timepoints))
  names(all_errors) <- timepoints

  result_models <- list()
  result <- list()
  for (k in timepoints) {
    excode_model <- excode_model_init

    if (verbose) {
      cat("Fitting model at position ", k, "\n", sep = "")
    }

    excode_model@transitions <- transMat_init
    excode_model@initial_prob <- initProb_init

    excode_model_fit <- NULL
    modelData <- NULL
    curr_error <- tryCatch(
      {
        modelData <- prepareData(surv_ts, excode_model, k,
          id = "ts1", time_units_back = time_units_back,
          past_weeks_not_included_state = past_timepoints_not_incuded,
          past_weeks_not_included_init = past_timepoints_not_incuded
        )

        if (learning_type %in% c("unsupervised", "semisupervised")) {
          excode_model_fit <- fitUnsupervised(excode_model, modelData, transMat_init,
            learning_type = learning_type, maxIter,
            verbose, time_units_back = time_units_back
          )
        } else if (learning_type == "supervised") {
          excode_model_fit <- fitSupervised(excode_model, modelData)
        } else {
          stop("learning_type must be one of: unsupervised, semisupervised or supervised.\n")
        }


        NA
      },
      error = function(e) {
        warning("Error fitting model at position ", k, ": ", as.character(e), ".")
        paste("Error fitting model at position ", k, ": ", as.character(e), ".", sep = "")
      }
    )


    if (!is.na(curr_error)) {
      ts_len <- 1
      if (is.null(excode_model_fit)) {
        excode_model_fit <- list(hmm = excode_model_init)
      }
      if (!is.null(modelData)) {
        ts_len <- length(which(excode_model_fit$model$state == 0))
      }

      nStates <- excode_model_fit$hmm@nStates
      excode_model_fit$hmm@posterior <- matrix(rep(NA, ts_len * nStates), ncol = nStates)
      excode_model_fit$hmm@alpha <- matrix(rep(NA, 2 * ts_len * nStates), ncol = nStates)

      excode_model_fit$hmm@converged <- rep(FALSE, ts_len)
      excode_model_fit$hmm@niter <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@LogLik <- rep(as.numeric(NA), ts_len)

      excode_model_fit$hmm@id <- rep(as.character(NA), ts_len)

      excode_model_fit$hmm@BIC <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@AIC <- rep(as.numeric(NA), ts_len)

      excode_model_fit$hmm@timepoint_fit <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@timepoint <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@date <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@observed <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@emission@mu0 <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@emission@mu1 <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@emission@mu <- matrix(rep(as.numeric(NA), ts_len * nStates), ncol = nStates)

      excode_model_fit$hmm@population <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@error <- rep(curr_error, ts_len)
      excode_model_fit$hmm@method <- rep(as.character(NA), ts_len)

      nb_size_NA <- data.frame(nb_size = rep(as.numeric(NA), ts_len))
      excode_model_fit$hmm@emission@distribution <-
        set_nb_size(
          excode_model_fit$hmm@emission@distribution,
          nb_size_NA, 1:nrow(nb_size_NA)
        )


      excode_model_fit$hmm@pval <- rep(as.numeric(NA), ts_len)
    } else {
      nStates <- length(unique(excode_model_fit$model$state))
      index_0 <- which(excode_model_fit$model$state == 0)
      index_1 <- which(excode_model_fit$model$state == 1)
      excode_model_fit$hmm@timepoint_fit <- rep(k, length(index_0))
      ts_len <- which(excode_model_fit$model$rtime == k &
        excode_model_fit$model$state == 0)
      excode_model_fit$hmm@posterior <- excode_model_fit$hmm_expectation$gamma

      y <- excode_model_fit$model$response[1:(nrow(excode_model_fit$model) / excode_model_fit$hmm@nStates)]
      rtime <- excode_model_fit$model$rtime[1:(nrow(excode_model_fit$model) / excode_model_fit$hmm@nStates)]
      transMat <- excode_model_fit$hmm@transitions
      initProb <- excode_model_fit$hmm@initial_prob
      states <- excode_model_fit$model$known_state[1:(nrow(modelData) / excode_model_fit$hmm@nStates)]
      emissionProb <- calcEmissionProb(
        excode_model_fit$hmm@emission@distribution,
        excode_model_fit$model
      )

      alpha_beta <- getAlphaBeta(y, transMat, initProb, emissionProb, states)
      rownames(alpha_beta$alpha) <- excode_model_fit$model$id[index_0]
      alpha_beta$alpha <- cbind(alpha_beta$alpha, rtime)
      alpha_beta$alpha <- cbind(alpha_beta$alpha, rep(k, nrow(alpha_beta$alpha)))
      # colnames(alpha_beta$alpha) <- c("alpha_normal", "alpha_excess", "rtime", "timepoint")
      excode_model_fit$hmm@alpha <- do.call("rbind", lapply(ts_len, function(x) {
        alpha_beta$alpha[(x - 1):x, ]
      }))


      excode_model_fit$hmm@converged <- rep(
        excode_model_fit$niter <= maxIter,
        length(index_0)
      )
      excode_model_fit$hmm@niter <- rep(excode_model_fit$niter, length(index_0))
      excode_model_fit$hmm@LogLik <- rep(excode_model_fit$loglik, length(index_0))


      excode_model_fit$hmm@id <- excode_model_fit$model$id[index_0]

      npar <- length(excode_model_fit$hmm@emission@excode_formula@coefficients) +
        getSharedParN(excode_model_fit$hmm@emission@distribution) *
          length(unique(excode_model_fit$hmm@id)) + 6
      excode_model_fit$hmm@BIC <- rep(log(nrow(excode_model_fit$hmm@posterior)) * npar -
        2 * excode_model_fit$hmm@LogLik[1], length(index_0))
      excode_model_fit$hmm@AIC <- rep(
        -2 * npar - 2 * excode_model_fit$hmm@LogLik[1],
        length(index_0)
      )

      excode_model_fit$hmm@timepoint <- excode_model_fit$model$rtime[index_0]
      excode_model_fit$hmm@date <- excode_model_fit$model$date[index_0]
      excode_model_fit$hmm@observed <- excode_model_fit$model$response[index_0]
      excode_model_fit$hmm@emission@mu0 <- excode_model_fit$model$mu[index_0]
      excode_model_fit$hmm@emission@mu1 <- excode_model_fit$model$mu[index_1]
      excode_model_fit$hmm@emission@mu <- matrix(excode_model_fit$model$mu, ncol = nStates)
      excode_model_fit$hmm@population <- excode_model_fit$model$population[index_0]
      excode_model_fit$hmm@error <- rep(as.character(curr_error), length(index_0))


      if (excode_model_fit$hmm@nStates == 2 &
        any(excode_model_fit$hmm@emission@mu0 > excode_model_fit$hmm@emission@mu1)) {
        excode_model_fit$hmm@emission@mu0 <- excode_model_fit$model$mu[index_1]
        excode_model_fit$hmm@emission@mu1 <- excode_model_fit$model$mu[index_0]
        excode_model_fit$hmm@emission@mu <- excode_model_fit$hmm@emission@mu[, 2:1]
        excode_model_fit$hmm@alpha <- excode_model_fit$hmm@alpha[2:1, ]
        excode_model_fit$hmm@posterior <- excode_model_fit$hmm@posterior[, 2:1]
        excode_model_fit$hmm@transitions <- excode_model_fit$hmm@transitions[2:1, 2:1]
        excode_model_fit$hmm@initial_prob <- excode_model_fit$hmm@initial_prob[2:1]
      }

      excode_model_fit$hmm@emission@distribution <-
        set_nb_size(
          excode_model_fit$hmm@emission@distribution,
          excode_model_fit$model, index_0
        )


      excode_model_fit$hmm@pval <- calculatePvalue(
        excode_model_fit$hmm@emission@distribution,
        excode_model_fit$hmm
      )
      excode_model_fit$hmm@method <- rep(learning_type, length(index_0))
    }


    if (!return_full_model) {
      nStates <- excode_model_fit$hmm@nStates
      ts_len <- which(excode_model_fit$hmm@timepoint == k &
        excode_model_fit$hmm@timepoint_fit == k)

      excode_model_fit$hmm@timepoint_fit <- excode_model_fit$hmm@timepoint_fit[ts_len]
      excode_model_fit$hmm@posterior <- matrix(excode_model_fit$hmm@posterior[ts_len, ], ncol = nStates)
      colnames(excode_model_fit$hmm@posterior) <- paste0("posterior", 0:(nStates - 1))
      # excode_model_fit$hmm@alpha <-
      #  matrix(excode_model_fit$hmm@alpha[(ts_len-1):ts_len,], ncol=2)
      excode_model_fit$hmm@id <- excode_model_fit$hmm@id[ts_len]

      excode_model_fit$hmm@timepoint <- excode_model_fit$hmm@timepoint[ts_len]
      excode_model_fit$hmm@date <- excode_model_fit$hmm@date[ts_len]

      excode_model_fit$hmm@observed <- excode_model_fit$hmm@observed[ts_len]
      excode_model_fit$hmm@emission@mu0 <- excode_model_fit$hmm@emission@mu0[ts_len]
      excode_model_fit$hmm@emission@mu1 <- excode_model_fit$hmm@emission@mu1[ts_len]
      excode_model_fit$hmm@emission@mu <- matrix(excode_model_fit$hmm@emission@mu[ts_len, ],
        ncol = nStates
      )
      colnames(excode_model_fit$hmm@emission@mu) <- paste0("mu", 0:(nStates - 1))
      excode_model_fit$hmm@population <- excode_model_fit$hmm@population[ts_len]

      excode_model_fit$hmm@emission@distribution <-
        subset_nb_size(excode_model_fit$hmm@emission@distribution, ts_len)


      excode_model_fit$hmm@pval <- excode_model_fit$hmm@pval[ts_len]
      excode_model_fit$hmm@method <- excode_model_fit$hmm@method[ts_len]
      excode_model_fit$hmm@error <- excode_model_fit$hmm@error[ts_len]

      excode_model_fit$hmm@converged <- excode_model_fit$hmm@converged[ts_len]
      excode_model_fit$hmm@niter <- excode_model_fit$hmm@niter[ts_len]
      excode_model_fit$hmm@LogLik <- excode_model_fit$hmm@LogLik[ts_len]
      excode_model_fit$hmm@BIC <- excode_model_fit$hmm@BIC[ts_len]
      excode_model_fit$hmm@AIC <- excode_model_fit$hmm@AIC[ts_len]
    }

    result_models[[as.character(k)]] <- excode_model_fit$hmm
  }

  merge_results(result_models)
}




merge_results <- function(excode_model_fit_list) {
  excode_model_fit <- excode_model_fit_list[[1]]

  if (length(excode_model_fit_list) > 1) {
    for (i in 2:length(excode_model_fit_list)) {
      excode_model_fit@timepoint_fit <- c(
        excode_model_fit@timepoint_fit,
        excode_model_fit_list[[i]]@timepoint_fit
      )
      excode_model_fit@posterior <- rbind(
        excode_model_fit@posterior,
        excode_model_fit_list[[i]]@posterior
      )
      excode_model_fit@alpha <- rbind(
        excode_model_fit@alpha,
        excode_model_fit_list[[i]]@alpha
      )

      excode_model_fit@converged <- c(
        excode_model_fit@converged,
        excode_model_fit_list[[i]]@converged
      )
      excode_model_fit@niter <- c(
        excode_model_fit@niter,
        excode_model_fit_list[[i]]@niter
      )
      excode_model_fit@LogLik <- c(
        excode_model_fit@LogLik,
        excode_model_fit_list[[i]]@LogLik
      )

      excode_model_fit@id <- c(
        excode_model_fit@id,
        excode_model_fit_list[[i]]@id
      )
      excode_model_fit@BIC <- c(
        excode_model_fit@BIC,
        excode_model_fit_list[[i]]@BIC
      )
      excode_model_fit@AIC <- c(
        excode_model_fit@AIC,
        excode_model_fit_list[[i]]@AIC
      )

      excode_model_fit@date <- c(
        excode_model_fit@date,
        excode_model_fit_list[[i]]@date
      )
      excode_model_fit@timepoint <- c(
        excode_model_fit@timepoint,
        excode_model_fit_list[[i]]@timepoint
      )
      excode_model_fit@observed <- c(
        excode_model_fit@observed,
        excode_model_fit_list[[i]]@observed
      )

      excode_model_fit@emission@mu0 <- c(
        excode_model_fit@emission@mu0,
        excode_model_fit_list[[i]]@emission@mu0
      )
      excode_model_fit@emission@mu1 <- c(
        excode_model_fit@emission@mu1,
        excode_model_fit_list[[i]]@emission@mu1
      )
      excode_model_fit@emission@mu <- rbind(
        excode_model_fit@emission@mu,
        excode_model_fit_list[[i]]@emission@mu
      )

      excode_model_fit@population <- c(
        excode_model_fit@population,
        excode_model_fit_list[[i]]@population
      )
      excode_model_fit@error <- c(
        excode_model_fit@error,
        excode_model_fit_list[[i]]@error
      )

      excode_model_fit@pval <- c(
        excode_model_fit@pval,
        excode_model_fit_list[[i]]@pval
      )
      excode_model_fit@method <- c(
        excode_model_fit@method,
        excode_model_fit_list[[i]]@method
      )

      excode_model_fit@emission@distribution <-
        merge_excode_family_result(
          excode_model_fit@emission@distribution,
          excode_model_fit_list[[i]]@emission@distribution
        )
    }
  }

  excode_model_fit
}




#' @title Detect Excess Counts in Epidemiological Time Series
#'
#' @description
#' This function fits an `excodeModel` to epidemiological surveillance data
#' and returns the estimated parameters of an `excodeFamily`. It can handle
#' single or multiple time series.
#'
#' @param surv_ts A surveillance time series, which can be a \code{data.frame},
#' an \code{sts} object, or a list of \code{sts} objects.
#' @param excode_model An object of class \code{\linkS4class{excodeModel}} which specifies the model parameters and structure
#' @param timepoints An \code{integer} or a sequence of integers specifying
#' the time points for which excess count detection should be performed.
#' @param learning_type  A \code{character} string indicating the type of learning.
#' Must be one of \code{c("unsupervised", "semisupervised", "supervised")}.
#' @param maxIter An \code{integer} specifying the maximum number of iterations
#' for the Expectation-Maximization (EM) algorithm. Defaults to \code{100}.
#' @param verbose A \code{logical} indicating whether progress should be printed
#' during execution. Defaults to \code{FALSE}.
#' @param return_full_model A \code{logical} indicating whether to return the full
#' model output. If \code{FALSE} (default), only results for the specified
#' \code{timepoints} are returned. If \code{TRUE}, the complete time series used
#' for model fitting is returned.
#' @param past_timepoints_not_incuded An \code{integer} specifying the number of past
#' time points to exclude from initialization. Defaults to \code{26}.
#' @param time_units_back Number of years to be used for model fitting.
#' @param timepoints_per_unit Number of time points within the considered time unit (e.g. 52 for weekly observations in a year).
#'
#' @returns Returns a fitted \code{\linkS4class{excodeModel}} object.
#'
#' @seealso \code{\linkS4class{excodeModel}}
#' @rdname run_excode
#'
#' @examples
#'
#' # Example 1: Run with harmonic Poisson model using data.frame
#' excode_family_pois <- excodeFamily("Poisson")
#' excode_formula_har <- excodeFormula("Harmonic")
#' excode_har_pois <- excodeModel(
#'   excode_family_pois,
#'   excode_formula_har
#' )
#' data(shadar_df)
#' result_shadar_har <- run_excode(shadar_df, excode_har_pois, 209:295)
#'
#'
#' # Example 2: Run with sts object (from 'surveillance' package)
#' \dontrun{
#' library(surveillance)
#' data(stsNewport)
#' result_newport_har <- run_excode(stsNewport, excode_har_pois, 209:295)
#' }
#'
#' # TODO need an example different usage of time_units_back and
#' # timepoints_per_unit which is simpler than the one in vignette
#'
#' @export
setGeneric("run_excode", function(surv_ts, excode_model, timepoints, learning_type = "unsupervised", maxIter = 100, verbose = FALSE, return_full_model = FALSE, past_timepoints_not_incuded = 26, time_units_back = 5, timepoints_per_unit = 52) standardGeneric("run_excode"))


#' Method for sts input
#' Direct application of run_excode
#' @rdname run_excode
setMethod("run_excode",
  signature = c("sts", "excodeModel"),
  function(surv_ts, excode_model, timepoints, learning_type = "unsupervised", maxIter = 100, verbose = FALSE, return_full_model = FALSE, past_timepoints_not_incuded = 26, time_units_back = 5, timepoints_per_unit = 52) {
    run_excode_internal(surv_ts, excode_model, timepoints, learning_type, maxIter, verbose, return_full_model, past_timepoints_not_incuded, time_units_back)
  }
)

#' Method for list of sts input
#' Sorting the sts objects by name
#' @rdname run_excode
setMethod("run_excode",
  signature = c("list", "excodeModel"),
  function(surv_ts, excode_model, timepoints, learning_type = "unsupervised", maxIter = 100, verbose = FALSE, return_full_model = FALSE, past_timepoints_not_incuded = 26, time_units_back = 5, timepoints_per_unit = 52) {
    run_excode_internal(surv_ts[sort(names(surv_ts))], excode_model, timepoints, learning_type, maxIter, verbose, return_full_model, past_timepoints_not_incuded, time_units_back)
  }
)
#' Method for data.frame input
#' Transformation of data.frame to sts object
#' @rdname run_excode
setMethod("run_excode",
  signature = c("data.frame", "excodeModel"),
  function(surv_ts, excode_model, timepoints, learning_type = "unsupervised", maxIter = 100, verbose = FALSE, return_full_model = FALSE, past_timepoints_not_incuded = 26, time_units_back = 5, timepoints_per_unit = 52) {
    if (!any(names(surv_ts) == "state")) {
      if (learning_type %in% c("semisupervised", "supervised")) {
        stop("Variable \'state\' must be provided with method \'", learning_type, "\'\n")
      }
      surv_ts$state <- NA
    }
    surv_ts_list <- split(surv_ts, surv_ts$id)
    for (i in 1:length(surv_ts_list)) {
      offset <- matrix(rep(1, nrow(surv_ts_list[[i]])), ncol = 1)
      if ("offset" %in% names(surv_ts_list[[i]])) {
        offset <- matrix(surv_ts_list[[i]]$offset, ncol = 1)
      }
      surv_ts_list[[i]] <- sts(surv_ts_list[[i]]$observed,
        frequency = timepoints_per_unit,
        epoch = surv_ts_list[[i]]$date,
        state = surv_ts_list[[i]]$state,
        population = offset
      )
    }
    run_excode_internal(
      surv_ts_list[sort(names(surv_ts_list))],
      excode_model, timepoints, learning_type,
      maxIter, verbose, return_full_model,
      past_timepoints_not_incuded, time_units_back
    )
  }
)
