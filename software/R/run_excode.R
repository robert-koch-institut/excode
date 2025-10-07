#' @title Detect Excess Counts in Epidemiological Time Series
#'
#' @description
#' Fit an \code{\linkS4class{excodeModel}} to epidemiological surveillance data
#' and return a fitted model (with posterior probabilities, fitted means, and
#' related diagnostics) for detecting excess counts. Supports single or multiple
#' time series provided as a \code{data.frame} only.
#'
#' @details
#' The function currently accepts a single \code{data.frame} as the surveillance
#' time series input. The \code{data.frame} must contain:
#' \itemize{
#'   \item \strong{date}: aggregation date;
#'   \item \strong{observed}: observed case counts;
#'   \item \strong{id}: identifier of the time series (use a single id for one series);
#'   \item \strong{offset} (optional): e.g., population at risk;
#'   \item \strong{state} (optional): periods with special events (e.g., 0 = normal,
#'         1 = excess, \code{NA} = unknown).
#' }
#'
#' It is no longer possible to pass \code{sts} objects (from the
#' \pkg{surveillance} package) or lists of such objects. These will result in an
#' error. A matrix-like collection of \code{sts} objects (e.g.,
#' \code{surveillance::measlesDE}) is also not supported.
#'
#' If the model's emission uses a \code{MultiState} \code{excodeFormula} that
#' already contains a \code{surv_ts} slot, \code{run_excode()} will use that data
#' and (if not supplied) set \code{timepoints} accordingly.
#'
#' @param excode_model An object of class \code{\linkS4class{excodeModel}}
#' specifying the model family and formula as well as initial parameters/structure.
#' @param surv_ts A \code{data.frame} as specified under \emph{Details}.
#' @param timepoints An \code{integer} or integer vector specifying the time points
#' at which excess-count detection should be performed. If \code{NULL} and the model
#' embeds a \code{surv_ts} via a \code{MultiState} formula, an appropriate default is chosen.
#' @param maxIter \code{integer}. Maximum number of EM iterations. Defaults to \code{100}.
#' @param verbose \code{logical}. Print progress during execution. Defaults to \code{FALSE}.
#' @param return_full_model \code{logical}. If \code{FALSE} (default), return only results
#' for the specified \code{timepoints}. If \code{TRUE}, return the complete time series used
#' for fitting.
#' to exclude from initialization. Defaults to \code{26}. For daily data this excludes the past
#' \emph{x} days; for weekly data, the past \emph{x} weeks.
#' @param time_units_back \code{integer}. How many past time \emph{units} to include when fitting.
#' Default \code{5}. For example, if \code{timepoints_per_unit = 52} (weekly data, units = years)
#' in the \code{excodeFormula}, then \code{5} corresponds to the past 5 years; if
#' \code{timepoints_per_unit = 7} (daily data, units = weeks), then \code{5} corresponds to the
#' past 5 weeks.
#'
#' @returns A fitted \code{\linkS4class{excodeModel}} object containing, for the requested
#' time points (or the full series if \code{return_full_model = TRUE}), posterior state
#' probabilities, fitted means, p-values, Anscombe residuals, convergence information, and
#' information criteria.
#'
#' @seealso \code{\linkS4class{excodeModel}}, \code{\linkS4class{excodeFamily}}, \code{\linkS4class{excodeFormula}}
#'
#' @examples
#' # Create a Poisson harmonic model
#' excode_family_pois <- excodeFamily("Poisson")
#' excode_formula_har <- excodeFormula("Harmonic")
#' excode_har_pois <- excodeModel(excode_family_pois, excode_formula_har)
#'
#' # Example: data.frame as input
#' data(shadar_df)
#' result_shadar_har <- run_excode(excode_har_pois, shadar_df, 209:295)
#'
#' @rdname run_excode
#' @export
run_excode <- function(surv_ts, timepoints=NULL,
                       time_units_back = 5, distribution, 
                       states, time_trend="None", 
                       periodic_model="Harmonic", 
                       period_length=52, 
                       intercept=TRUE,
                       covariate_df=NULL,
                       weights_threshold = 2.58,
                       weights_threshold_baseline = 1, 
                       maxIter = 100,
                       verbose = FALSE, 
                       return_full_model = FALSE) {
  
  
  if(is.null(timepoints)) {
    stop("Specify timepoints for model fitting.")
  } else if(length(timepoints)>1 & return_full_model) {
    stop("`return_full_model==TRUE` not allowed with `length(timepoints)>0`.")
  }
  surv_ts$id <- "ts1"
  surv_ts_init <- surv_ts
  no_states <- states 
  
  all_errors <- rep("", length(timepoints))
  names(all_errors) <- timepoints

  result_models <- list()
  result <- list()
  for (k in timepoints) {
    excode_model <- init_excode(surv_ts_init, k, time_units_back, distribution, 
                                no_states, time_trend=time_trend, periodic_model, 
                               period_length=period_length, 
                               intercept=intercept,
                               covariate_df=covariate_df,
                               weights_threshold = weights_threshold,
                               weights_threshold_baseline = weights_threshold_baseline) 
    
    v <- validate_run_excode_inputs(excode_model, nrow(excode_model@emission@excode_formula@surv_ts),
                                    maxIter, verbose, 
                                    return_full_model, time_units_back)
    surv_ts  <- v$surv_ts
    timepoints <- v$timepoints
    maxIter  <- v$maxIter
    verbose  <- v$verbose
    return_full_model <- v$return_full_model
    time_units_back <- v$time_units_back
    
    transMat_init <- excode_model@transitions
    initProb_init <- excode_model@initial_prob
    
   
    if (verbose) {
      cat("Fitting model at position ", k, "\n", sep = "")
    }

    excode_model_fit <- NULL
    modelData <- NULL
    curr_error <- tryCatch(
      {
        modelData <- prepareData(surv_ts, excode_model, nrow(surv_ts),
          id = "ts1", time_units_back = time_units_back
        )

          excode_model_fit <- fitUnsupervised(excode_model, modelData, transMat_init,
             maxIter, verbose, time_units_back = time_units_back)
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
        excode_model_fit <- list(hmm = excode_model)
      }
      ts_len <- nrow(excode_model_fit$hmm@emission@excode_formula@surv_ts)
      
      nStates <- excode_model_fit$hmm@nStates
      posterior <- matrix(rep(NA, ts_len * nStates), ncol = nStates)
      colnames(posterior) <- paste0("posterior", 0:(nStates-1))
      excode_model_fit$hmm@posterior <- posterior
      excode_model_fit$hmm@alpha <- matrix(rep(NA, 2 * ts_len * nStates), ncol = nStates)

      excode_model_fit$hmm@converged <- rep(FALSE, ts_len)
      excode_model_fit$hmm@niter <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@LogLik <- rep(as.numeric(NA), ts_len)

      excode_model_fit$hmm@id <- rep(as.character(NA), ts_len)

      excode_model_fit$hmm@BIC <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@AIC <- rep(as.numeric(NA), ts_len)

      excode_model_fit$hmm@timepoint_fit <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@timepoint <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@date <- rep(as.Date(NA), ts_len)
      excode_model_fit$hmm@observed <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@emission@mu0 <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@emission@mu1 <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@emission@mu <- matrix(rep(as.numeric(NA), ts_len * nStates), ncol = nStates)

      excode_model_fit$hmm@population <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@error <- rep(curr_error, ts_len)
      
      nb_size_NA <- data.frame(nb_size = rep(as.numeric(NA), ts_len))
      excode_model_fit$hmm@emission@distribution <-
        set_nb_size(
          excode_model_fit$hmm@emission@distribution,
          nb_size_NA, 1:nrow(nb_size_NA)
        )


      excode_model_fit$hmm@pval <- rep(as.numeric(NA), ts_len)
      excode_model_fit$hmm@anscombe_residual <- rep(as.numeric(NA), ts_len)
    } else {
      nStates <- length(unique(excode_model_fit$model$state))
      index_0 <- which(excode_model_fit$model$state == 0)
      index_1 <- which(excode_model_fit$model$state == 1)
      excode_model_fit$hmm@timepoint_fit <- rep(k, length(index_0))
      ts_len <- which(excode_model_fit$model$rtime == nrow(surv_ts) &
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

      excode_model_fit$hmm@timepoint <- (k-length(index_0)+1):k # excode_model_fit$model$rtime[index_0]
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
      
      excode_model_fit$hmm@anscombe_residual <- anscombe_residuals(
        excode_model_fit$hmm@emission, 
        excode_model_fit$model
      )
    }


    if (!return_full_model) {
      nStates <- excode_model_fit$hmm@nStates
      ts_len <- which(excode_model_fit$hmm@timepoint == nrow(surv_ts) &
        excode_model_fit$hmm@timepoint_fit == k)
      excode_model_fit$hmm@timepoint_fit <- excode_model_fit$hmm@timepoint_fit[ts_len]
      excode_model_fit$hmm@posterior <- matrix(excode_model_fit$hmm@posterior[ts_len, ], ncol = nStates)
      colnames(excode_model_fit$hmm@posterior) <- paste0("posterior", 0:(nStates - 1))
      # excode_model_fit$hmm@alpha <-
      #  matrix(excode_model_fit$hmm@alpha[(ts_len-1):ts_len,], ncol=2)
      excode_model_fit$hmm@id <- excode_model_fit$hmm@id[ts_len]

      excode_model_fit$hmm@timepoint <- k#excode_model_fit$hmm@timepoint[ts_len]
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
      excode_model_fit$hmm@anscombe_residual <- excode_model_fit$hmm@anscombe_residual[ts_len]
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
      
      excode_model_fit@anscombe_residual <- c(
        excode_model_fit@anscombe_residual,
        excode_model_fit_list[[i]]@anscombe_residual
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




validate_run_excode_inputs <- function(excode_model,
                                       timepoints,
                                       maxIter,
                                       verbose,
                                       return_full_model,
                                       time_units_back) {
  
  surv_ts <- excode_model@emission@excode_formula@surv_ts
  
  # 1) excode_model -----------------------------------------------------------
  if (!methods::is(excode_model, "excodeModel")) {
    stop("`excode_model` must be an S4 object of class 'excodeModel'.")
  }
  required_slots <- c("emission", "transitions", "initial_prob")
  miss_slots <- setdiff(required_slots, methods::slotNames(excode_model))
  if (length(miss_slots)) {
    stop("`excode_model` is missing required slot(s): ",
         paste(miss_slots, collapse = ", "))
  }
  
  # 2) surv_ts: only data.frame allowed (no sts or lists of sts) --------------
  if (!is.null(surv_ts)) {
    if (methods::is(surv_ts, "sts")) {
      stop("`surv_ts` may not be an 'sts' object. Provide a data.frame instead.")
    }
    if (is.list(surv_ts) && any(vapply(surv_ts, methods::is, logical(1), "sts"))) {
      stop("Lists or matrices of 'sts' objects are not supported. Provide a data.frame.")
    }
    if (!is.data.frame(surv_ts)) {
      stop("`surv_ts` must be a data.frame (or NULL if provided via MultiState formula).")
    }
  }
  
  # 3) timepoints --------------------------------------------------------------
  if (!is.null(timepoints)) {
    if (!is.numeric(timepoints)) stop("`timepoints` must be numeric (integer-like).")
    if (anyNA(timepoints)) stop("`timepoints` may not contain NA.")
    if (any(timepoints <= 0)) stop("`timepoints` must be positive (1-based indices).")
    if (!all(timepoints == as.integer(timepoints))) {
      warning("Coercing `timepoints` to integers.")
      timepoints <- as.integer(round(timepoints))
    } else {
      timepoints <- as.integer(timepoints)
    }
    timepoints <- sort(unique(timepoints))
  }
  
  # 4) scalars / types --------------------------------------------------------
  chk_scalar_logical <- function(x, nm) {
    if (!is.logical(x) || length(x) != 1 || is.na(x)) {
      stop("`", nm, "` must be a single non-NA logical (TRUE/FALSE).")
    }
  }
  chk_scalar_int <- function(x, nm, min_val = -Inf, allow_na = FALSE) {
    if (length(x) != 1 || (!allow_na && is.na(x))) {
      stop("`", nm, "` must be a single value.")
    }
    if (!is.numeric(x)) stop("`", nm, "` must be numeric (integer-like).")
    if (x != as.integer(x)) {
      warning("Coercing `", nm, "` to integer.")
      x <<- as.integer(round(x))
    } else {
      x <<- as.integer(x)
    }
    if (!is.infinite(min_val) && x < min_val) {
      stop("`", nm, "` must be >= ", min_val, ".")
    }
    return(x)
  }
  
  maxIter                      <- chk_scalar_int(maxIter, "maxIter", min_val = 1)
  time_units_back              <- chk_scalar_int(time_units_back, "time_units_back",
                                                 min_val = 1)
  chk_scalar_logical(verbose, "verbose")
  chk_scalar_logical(return_full_model, "return_full_model")
  

  # 6) surv_ts column checks & normalization ----------------------------------
  normalize_surv_ts <- function(df) {
    req_cols <- c("date", "observed")
    miss <- setdiff(req_cols, names(df))
    if (length(miss)) {
      stop("`surv_ts` is missing required column(s): ", paste(miss, collapse = ", "))
    }
    
    # date
    if (!inherits(df$date, "Date")) {
      if (inherits(df$date, "POSIXt")) {
        df$date <- as.Date(df$date)
      } else {
        suppressWarnings({
          try_date <- as.Date(df$date)
        })
        if (any(is.na(try_date))) {
          stop("`surv_ts$date` must be Date-like. Failed to convert some entries.")
        }
        df$date <- try_date
      }
    }
    
    # observed
    if (!is.numeric(df$observed)) stop("`surv_ts$observed` must be numeric.")
    if (anyNA(df$observed)) stop("`surv_ts$observed` contains NA.")
    if (any(df$observed < 0)) stop("`surv_ts$observed` must be >= 0.")
    if (any(df$observed != as.integer(df$observed))) {
      warning("`surv_ts$observed` is not integer; rounding to nearest integer.")
      df$observed <- as.integer(round(df$observed))
    }
    
    
    # offset
    if (!("offset" %in% names(df))) {
      message("No `offset` column; assuming offset = 1 for all rows.")
      df$offset <- 1
    }
    if (!is.numeric(df$offset) || anyNA(df$offset) || any(df$offset < 1)) {
      stop("`surv_ts$offset` must be >1 without NA.")
    }
    
    # ---- STATE: strict rule: numeric and only 0 or NA -----------------------
    if (!("state" %in% names(df))) {
      message("No `state` column; creating `state = NA`.")
      df$state <- NA_real_
    }
    # must be numeric, no coercions performed
    if (!is.numeric(df$state)) {
      stop("`surv_ts$state` must be numeric and contain only 0 or NA.")
    }
    # allow NA; all non-NA must be exactly 0
    bad_state <- !is.na(df$state) & df$state != 0
    if (any(bad_state)) {
      stop("`surv_ts$state` must contain only 0 or NA (found: ",
           paste(unique(df$state[bad_state]), collapse = ", "), ").")
    }
    # store as integer for consistency
    df$state <- as.integer(ifelse(is.na(df$state), NA_integer_, 0L))
    # -------------------------------------------------------------------------
    
    # duplicates & order
    if (anyDuplicated(df[c("date")])) {
      dup <- df[duplicated(df[c("date")]) | duplicated(df[c("date")], fromLast = TRUE),
                c("date")]
      stop("`surv_ts` has duplicate (id, date) rows. Example: ",
           paste(utils::head(paste(dup$id, dup$date), 3), collapse = "; "), " …")
    }
    o <- order(df$id, df$date)
    if (!identical(o, seq_len(NROW(df)))) {
      warning("Sorting `surv_ts` by id, then date.")
      df <- df[o, , drop = FALSE]
    }
    
    # extra columns
    extra <- setdiff(names(df), c("date","observed","offset","state"))
    if (length(extra)) {
      message("Additional columns in `surv_ts` will be passed through/ignored as appropriate: ",
              paste(extra, collapse = ", "))
    }
    
    df
  }
  
  if (!is.null(surv_ts)) {
    surv_ts <- normalize_surv_ts(surv_ts)
  }
  
  # 7) timepoints vs data length ----------------------------------------------
  if (!is.null(timepoints) && !is.null(surv_ts)) {
    n <- NROW(surv_ts)
    if (any(timepoints > n)) {
      stop("Some `timepoints` (", paste(timepoints[timepoints > n], collapse = ", "),
           ") exceed the number of rows in `surv_ts` (", n, ").")
    }
  }
  
  # 8) minimal history feasibility checks -------------------------------------
  tpu <- NA_integer_
  try({
    tpu <- excode_model@emission@excode_formula@timepoints_per_unit
  }, silent = TRUE)
  if (!is.null(surv_ts) && length(tpu) == 1 && is.finite(tpu) && !is.na(tpu) && tpu > 0) {
    required_history <- time_units_back * as.integer(tpu)
    required_prefix  <- required_history
    if (!is.null(timepoints)) {
      bad <- timepoints <= required_prefix
      if (any(bad)) {
        warning("Some `timepoints` (", paste(timepoints[bad], collapse = ", "),
                ") may have insufficient history for fitting (need > ",
                required_prefix, " observations). Proceeding with available data.")
      }
    }
  }
  
  # 9) Return normalized values ----------------------------------------------
  list(
    surv_ts  = surv_ts,
    timepoints = timepoints,
    maxIter = maxIter,
    verbose = verbose,
    return_full_model = return_full_model,
    time_units_back = time_units_back
  )
}



