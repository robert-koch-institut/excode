#' This class is a generic container for the model used for outbreak detection
#'
#' @slot nStates Number of states in the hidden Markov model.
#' @slot emission Emission of hidden Markov model. Includes a excodeFamily and a excodeFormula object.
#' @slot transitions Transition proabilities of hidden Markov model.
#' @slot transitions_prior Indicates whether a prior for transition probabilities should be used.
#' @slot prior_weights Dirichlet prior weights for transition probabilities. This cannot be changed by the user.
#' @slot loglik_transitions Prior log likelihood of transition probabilities. Only relevant for model fitting.
#' @slot initial_prob Initial state probabilities of the hidden Markov model.
#' @slot setBckgState Indicates whether a background state should be computed for model fitting. Background states are set to 0s and time points where the Anscombe residuals of the initial model are < 1.
#' @slot id Name of the time series
#' @slot converged Logical, indicating whether the EM-algorithm converged.
#' @slot niter Maximum number of iterations for model fitting.
#' @slot LogLik Log lokelihood of the fitted model.
#' @slot AIC Akaike information criterion  of the fitted model.
#' @slot BIC Bayesian information criterion  of the fitted model.
#' @slot posterior Matrix of posterior probabilities for each time point (rows) and state (columns).
#' @slot alpha Alpha (Forward) probabilities of the HMM.
#' @slot pval P-value for testing wheter timepoint is in normal or excess state.
#' @slot timepoint_fit Numeric time point in timeseries used for fitting.
#' @slot timepoint Numeric time point in timeseries used for fitting.
#' @slot date Date for each time point.
#' @slot observed Vector of observed number of cases.
#' @slot population Population size, default: 1.
#' @slot error Character string of error message (if an error occured).
#'
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
    id = "character",
    converged = "logical",
    niter = "numeric",
    LogLik = "numeric",
    AIC = "numeric",
    BIC = "numeric",
    posterior = "matrix",
    alpha = "matrix",
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
#' @examples
#' # Initialisation of a mean model without timetrend with Poisson emission
#'
#' excode_formula_mean <- excodeFormula("Mean", timeTrend = FALSE)
#' excode_family_pois <- excodeFamily("Poisson")
#' excodeModel(excode_family_pois, excode_formula_mean)
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


#' @title Summary of an excodeModel.
#'
#' @description
#' Provides a summary of the \code{excodeModel} object, containing the results for those time points where excess count detection was performed.
#' By default, the summary includes model-based expectations, posterior probabilities, p-values, and their corresponding signal thresholds,
#' computed according to the model specification.
#'
#' @param object An object of class \code{excodeModel}. The fitted model to summarize.
#' @param pars Character vector. Specifies which parameters or variables to extract and include in the summary.
#'        Typical entries include: \code{"posterior"}, \code{"pval"}, \code{"date"}, \code{"timepoint"},
#'        \code{"observed"}, \code{"emission"}, \code{"id"}, \code{"BIC"}, \code{"AIC"}. Default includes all.
#' @param prob_threshold Numeric. Posterior probability threshold used to calculate the upper bound for the expected number of cases under the posterior.
#'        This value controls the sensitivity of alarm detection: a lower \code{prob_threshold} results in a lower \code{posterior_ub}.
#'        The computed upper bound is returned in the \code{posterior_ub} column of the summary output.
#'        Observations with a posterior probability greater than or equal to this threshold are marked as excess count.
#'        Default is 0.5.
#' @param pval_threshold Numeric. p-value threshold used used to calculate the upper bound for the expected number of cases using quantiles based on the p-value.
#'        This threshold determines the sensitivity of detection: lower values result in more conservative signal classification.
#'        It is used to compute the \code{pval_ub} column in the summary output.
#'        Observations with a p-value less than or equal to this threshold are marked as excess counts.
#'        Default is 0.01.

#' @param maxiter Integer. Maximum number of iterations to use when estimating the posterior alarm threshold. Default is 1000.
#'
#' @return A data.frame summarizing the selected components of the \code{excodeModel}, including expected values, posterior,
#' p-value, and model fit metrics (such as AIC and BIC) depending on the values of \code{pars}.
#'
#' @seealso \code{\linkS4class{excodeModel}}, \code{\linkS4class{excodeFamily}}, \code{\linkS4class{excodeFormula}}
#'
#' @examples
#'
#' # Looking at summary of the results using a harmonic Poisson model on the shadar_df
#' \dontrun{
#' #' excode_family_pois <- excodeFamily("Poisson")
#' excode_formula_har <- excodeFormula("Harmonic")
#' excode_har_pois <- excodeModel(excode_family_pois, excode_formula_har)
#' # perform excess count detection for time points 209:295
#' result_shadar_har <- run_excode(shadar_df, excode_har_pois, 209:295)
#' # obtain the summary of the results for the time points 209:295
#' summary(result_shadar_har)
#' }

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






#' @title Summary of an excodeModel.
#'
#' @param x An excodeModel object.
#'
#' @seealso \code{\linkS4class{excodeModel}}, \code{\linkS4class{excodeFamily}}, \code{\linkS4class{excodeFormula}}
#'
#' @examples
#'
#' # TODO
#'
#' @export
setMethod(f = "plot", signature = "excodeModel", function(x) {
  excode:::plot_model(x)
})
