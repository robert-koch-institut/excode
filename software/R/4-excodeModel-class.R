#' @title excodeModel Class
#'
#' @description The `excodeModel` class is a generic container for the model used for outbreak detection based on a hidden Markov model framework.
#'
#' @slot nStates Number of states in the hidden Markov model.
#' @slot emission Emission component of the hidden Markov model, including an `excodeFamily` and an `excodeFormula` object.
#' @slot transitions Transition probabilities of the hidden Markov model.
#' @slot transitions_prior Logical indicating whether a prior for transition probabilities should be used.
#' @slot prior_weights Dirichlet prior weights for transition probabilities; cannot be changed by the user.
#' @slot loglik_transitions Prior log-likelihood of transition probabilities; only relevant for model fitting.
#' @slot initial_prob Initial state probabilities of the hidden Markov model.
#' @slot setBckgState Logical indicating whether a background state should be computed for model fitting; background states are set to 0 for time points where the Anscombe residuals of the initial model are < 1.
#' @slot converged Logical indicating whether the EM algorithm converged.
#' @slot niter Number of iterations used for model fitting.
#' @slot LogLik Log-likelihood of the fitted model.
#' @slot AIC Akaike information criterion of the fitted model.
#' @slot BIC Bayesian information criterion of the fitted model.
#' @slot posterior Matrix of posterior probabilities for each time point (rows) and state (columns).
#' @slot pval P-values for testing whether each time point is in a normal or excess state.
#' @slot zscore Standardized Anscombe residuals (z-scores) for each time point.
#' @slot timepoint_fit Numeric time point in the time series used for fitting.
#' @slot timepoint Numeric time point in the time series.
#' @slot date Date corresponding to each time point.
#' @slot observed Vector of observed counts or case numbers.
#' @slot population Population size; default is 1.
#' @slot error Character string containing an error message, if an error occurred during model fitting.
#'
#' @exportClass excodeModel
setClass("excodeModel",
  slots = c(
    nStates = "numeric",
    emission = "Emission",
    transitions = "matrix",
    transitions_prior = "logical",
    prior_weights = "matrix",
    loglik_transitions = "numeric",
    initial_prob = "numeric",
    setBckgState = "logical",
    converged = "logical",
    niter = "numeric",
    LogLik = "numeric",
    AIC = "numeric",
    BIC = "numeric",
    posterior = "matrix",
    pval = "numeric",
    zscore = "numeric",
    timepoint_fit = "numeric",
    timepoint = "numeric",
    date = "Date",
    observed = "numeric",
    population = "numeric",
    error = "character"
  )
)

#' Prints description of FarringtonNoufaily object
#'
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("excodeModel"), function(object) {
  cat(ifelse(all(is.na(object@converged)), "Inital ",
    "Fitted "
  ), "excodeModel\n", sep = "")
  cat(" ")
  print(object@emission@distribution)
  cat(" ")
  print(object@emission@excode_formula)
  cat("No. of timepoints: ", length(unique(object@timepoint_fit)),
    "\n",
    sep = ""
  )
})


#' @title Create a model for excess count detection
#'
#' @description
#' Creates an object of class \code{excodeModel} for modeling excess counts using
#' a Hidden Markov Model (HMM) with user-defined emission distributions
#' and formula for modeling observed counts.
#'
#' @param family An \code{\linkS4class{excodeFamily}} object defining the emission distribution
#'   (e.g., Poisson, Negative Binomial).
#' @param formula An \code{\linkS4class{excodeFormula}} object specifying the structure of the
#'   model (e.g., time trends, seasonality, ...).
#' @param initial_mu Initial estimates of the mean for 'MultiState' models.
#' @param transMat Inital transition probabilities.
#' @param initProb A numeric vector containing initial state probabilities (probabilities of of states at first time point) of the hidden Markov model.
#' @param transMat_prior Logical. Should a prior distribution be used for estimating transition
#'   probabilities? Default is \code{TRUE}.
#' @param setBckgState Logical. Should a background state be inferred for model fitting?
#'   Background states are initialized to 0 for time points with Anscombe residuals < 1
#'   from an initial model. Default is \code{TRUE}.
#'
#' @seealso
#' \code{\linkS4class{excodeModel}},
#' \code{\linkS4class{excodeFamily}},
#' \code{\linkS4class{excodeFormula}}
#'
#' @return An object of class \code{\linkS4class{excodeModel}}.
#'
#'
#' @export
excodeModel <- function(family, formula,
                        initial_mu = NULL,
                        transMat = NULL,
                        initProb = NULL,
                        transMat_prior = TRUE,
                        setBckgState = TRUE) {
  emission <- Emission(family, formula, initial_mu)
  nStates <- 2
  if (inherits(formula, what = "MultiState")) {
    nStates <- formula@nStates
  }

  if (nStates == 2) {
    if (is.null(transMat)) {
      transMat <- matrix(c(0.95, 0.4, 0.05, 0.6), ncol = 2)
    }
    if (is.null(initProb)) {
      initProb <- c(0.5, 0.5)
    }
  } else {
    if (is.null(transMat)) {
      transMat <- matrix(rep(initProb, length = nStates),
        ncol = nStates,
        nrow = nStates
      )
    }
    if (is.null(initProb)) {
      initProb <- rep(1 / nStates, length = nStates)
    }
    setBckgState <- FALSE
    transMat_prior <- TRUE
  }


  new("excodeModel",
    emission = emission,
    nStates = nStates,
    transitions = transMat,
    transitions_prior = transMat_prior,
    prior_weights = matrix(NA, nrow = nStates, ncol = nStates),
    loglik_transitions = as.numeric(0),
    initial_prob = initProb,
    setBckgState = setBckgState,
    date = as.Date(NA)
  )
}

updateTransInitProb <- function(excode_model, gamma, xsi, modelData) {
  transMat <- Reduce("+", xsi)
  gammaSum_ind <- which(modelData$timepoint != max(modelData$timepoint) & modelData$state == 0)
  # print(gammaSum_ind)
  gammaSum <- apply(gamma[gammaSum_ind, ], 2, sum)
  if (excode_model@transitions_prior) { # Update with prior weights
    log_dir <- 0
    for (i in 1:nrow(transMat)) {
      transMat[i, ] <- (excode_model@prior_weights[i, ] - 1 + transMat[i, ]) / (sum(excode_model@prior_weights[i, ]) - ncol(transMat) + gammaSum[i])
      log_dir <- log_dir + log_dirichlet(transMat[i, ], excode_model@prior_weights[i, ])
    }
  } else { # Update without prior
    for (i in 1:nrow(transMat)) {
      transMat[i, ] <- transMat[i, ] / gammaSum[i]
    }
  }

  initProb_ind <- which(modelData$timepoint == 0 & modelData$state == 0)
  initProb <- gamma[initProb_ind, ]
  if (length(initProb_ind) > 1) {
    initProb <- apply(initProb, 2, sum)
  }
  initProb <- initProb / sum(initProb)

  excode_model@transitions <- transMat
  excode_model@initial_prob <- as.vector(initProb)
  excode_model@loglik_transitions <- log_dir

  excode_model
}


#' @title Summary of an excodeModel
#'
#' @description Summarize a fitted \code{excodeModel} at time points where excess-count detection was performed, including expectations, posterior probabilities, p-values, and corresponding thresholds per the model specification.
#'
#' @param object An \code{excodeModel} to summarize.
#' @param pars Character vector of fields to include (e.g., \code{"posterior"}, \code{"pval"}, \code{"zscore"}, \code{"date"}, \code{"timepoint"}, \code{"observed"}, \code{"emission"}, \code{"BIC"}, \code{"AIC"}); defaults to all.
#' @param prob_threshold Numeric posterior probability threshold used to compute the posterior-based upper bound (\code{posterior_ub}) and flag excess counts (>= threshold); default 0.5.
#' @param pval_threshold Numeric p-value threshold used to compute the p-value-based upper bound (\code{pval_ub}) and flag excess counts (<= threshold); default 0.01.
#' @param anscombe_threshold Numeric threshold on the Anscombe z-score used for signal assessment; default 2.
#' @param maxiter Integer maximum iterations when estimating the posterior alarm threshold; default 1000.
#'
#' @return A \code{data.frame} summarizing selected components of the \code{excodeModel} (expected values, posterior, p-value, and fit metrics such as AIC/BIC) according to \code{pars}.
#'
#' @seealso \code{\linkS4class{excodeModel}}, \code{\linkS4class{excodeFamily}}, \code{\linkS4class{excodeFormula}}
#'
#' @examples
#'
#' data(shadar_df)
#' res_har_pois <- run_excode(surv_ts = shadar_df, timepoints = 295, distribution = "Poisson", 
#' states = 2, periodic_model = "Harmonic", time_trend = "Linear", set_baseline_state = TRUE)
#' summary(res_har_pois)
#'
#' @export
setMethod("summary",
  signature = c("excodeModel"),
  function(object, pars = c(
             "posterior", "pval", "zscore", "date", "timepoint",
             "observed", "emission", "BIC", "AIC"
           ),
           prob_threshold = 0.5,
           pval_threshold = 0.05,
           anscombe_threshold = 2,
           maxiter = 1000) {
    emission_df <- NULL
    posterior_df <- NULL
    remove_timpoint_fit <- FALSE
    if (!"timepoint_fit" %in% pars) {
      pars <- c(pars, "timepoint_fit")
      remove_timpoint_fit <- TRUE
    }
    pars_order <- pars
    if ("emission" %in% pars) {
      par_seq <- 1:length(pars)
      emission_pos <- which(pars == "emission")
      emission_df <- summary_emission(object@emission)
      emission_pos <- emission_pos + (0:(ncol(emission_df) - 1))
      par_seq <- par_seq[pars != "emission"]
      pars <- pars[pars != "emission"]
      while (any(emission_pos %in% par_seq)) {
        change_this <- which(par_seq >= min(emission_pos))
        par_seq[change_this] <- par_seq[change_this] + 1
      }
      pars_order <- c(pars, names(emission_df))[order(c(par_seq, emission_pos))]
    }
    if ("posterior" %in% pars) {
      par_seq <- 1:length(pars_order)
      posterior_pos <- which(pars_order == "posterior")
      posterior_df <- object@posterior
      posterior_pos <- posterior_pos + (0:(ncol(posterior_df) - 1))
      par_seq <- par_seq[pars_order != "posterior"]
      pars <- pars[pars != "posterior"]
      pars_order <- pars_order[pars_order != "posterior"]
      while (any(posterior_pos %in% par_seq)) {
        change_this <- which(par_seq >= min(posterior_pos))
        par_seq[change_this] <- par_seq[change_this] + 1
      }
      pars_order <- c(pars_order, colnames(posterior_df))[order(c(par_seq, posterior_pos))]
    }

    res_df <- as.data.frame(lapply(pars, function(x) slot(object, x)))
    names(res_df) <- pars

    if (!is.null(emission_df)) {
      res_df <- cbind(res_df, emission_df)
    }
    if (!is.null(posterior_df)) {
      res_df <- cbind(res_df, posterior_df)
    }


    res_df <- res_df[pars_order]


    if (remove_timpoint_fit) {
      pars_order <- pars_order[pars_order != "timepoint_fit"]
    }


    res_df[pars_order]
  }
)
