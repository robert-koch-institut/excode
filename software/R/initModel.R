init_glm_mutlistate <- function(hmm,
                                modelData,
                                setBckgState = TRUE) {
  nStates <- hmm@nStates
  ts_len <- nrow(modelData)
  dataGLM_expanded <- modelData[rep(1:nrow(modelData), nStates), ]
  subset_mu <- modelData$rtime
  dataGLM_expanded$mu <- as.vector(hmm@emission@mu[subset_mu, ])
  dataGLM_expanded$state <- sort(rep(0:(nStates - 1), nrow(modelData)))

  for (curr_state in 1:(nStates - 1)) {
    dataGLM_expanded[[paste0("state", curr_state)]] <- as.numeric(curr_state == dataGLM_expanded$state)
  }
  states_formula <- paste(paste0("state", 1:(nStates - 1)), collapse = " + ")
  hmm@emission@excode_formula@formula <- paste0(
    hmm@emission@excode_formula@formula_bckg,
    " + ", states_formula
  )

  dataGLM_expanded$id_state <- factor(paste0("s", dataGLM_expanded$state),
    levels = paste0("s", 0:(nStates - 1))
  )
  dataGLM_expanded$nb_size <- NA
  if (hmm@emission@distribution@name == "NegBinom") {
    nb_size <- hmm@emission@distribution@nb_size
    nb_size <- unlist(lapply(nb_size, function(x) rep(x, length(unique(modelData$rtime)))))
    nb_size <- rep(nb_size, nStates)
    dataGLM_expanded$nb_size <- nb_size
  }

  dataGLM_expanded$known_state <- dataGLM_expanded$true_state


  emissionProb <- calcEmissionProb(
    hmm@emission@distribution,
    dataGLM_expanded
  )

  gamma_xsi <- forwardBackward(
    hmm, dataGLM_expanded,
    emissionProb
  )

  omega <- as.vector(gamma_xsi$gamma)
  model_updated <- updateEmission(hmm@emission, dataGLM_expanded, omega)
  model_updated
}


#' @title Initialize a multi-state EXCODE model from a surveillance time series
#'
#' @description Builds an initial multi-state EXCODE model by fitting a baseline GLM to the background component, computing Anscombe residuals with Farrington-style down-weighting, clustering non-baseline residuals into `states - 1` groups to seed latent states, and estimating state-specific mean trajectories for `initial_mu`.
#'
#' @param surv_ts data.frame of surveillance counts (response and any covariates required by the EXCODE formula); if an `offset` column is present it is used when building the MultiState formula.
#' @param timepoint Integer index indicating the right end of the training window passed to `extractModelData()`.
#' @param time_units_back Positive integer giving how many past time units (relative to `timepoint`) to include for initialization.
#' @param distribution Character scalar passed to `excodeFamily()` defining the emission distribution, e.g., `"Poisson"` or `"NegBinom"` (NB size inferred from the baseline fit).
#' @param states Integer (>= 2) number of latent states; state `0` is treated as background.
#' @param time_trend Character time-trend setting forwarded to `excodeFormula(name = periodic_model, time_trend = ...)`, one of `c("Linear","Spline1","Spline2","Spline3","Spline4","None")`; default `"None"`.
#' @param periodic_model Character model name for the background formula, one of `c("Mean","FarringtonNoufaily","Harmonic","Custom")`.
#' @param period_length Integer mapped to `timepoints_per_unit` in `excodeFormula` (e.g., `52` weekly points per year, `7` daily with weekly cycle, `365` daily with annual cycle); default 52.
#' @param intercept Logical forwarded to `excodeFormula(name = "MultiState", intercept = ...)`; default `TRUE`.
#' @param covariate_df data.frame or `NULL`; optional covariates used when constructing `excodeFormula()`; defaults to `NULL`.
#' @param weights_threshold Numeric SD cutoff used by `surveillance::algo.farrington.assign.weights()` to down-weight large residuals in the baseline refit; default `2.58`.
#' @param weights_threshold_baseline Numeric Anscombe residual threshold for identifying non-baseline time points that are clustered to seed the `states - 1` non-baseline states; default `1`.
#' @param set_baseline_state Logical; if `TRUE`, computes/assigns the background state after initialization; default `FALSE`.
#'
#' @details Internally constructs a background `excodeFormula(name = periodic_model, ...)`, fits a quasi-Poisson GLM with Farrington-style weights, clusters residuals above `weights_threshold_baseline` via `kmeans()` to initialize non-background states, then fits a GLM with a `state` factor to obtain `initial_mu`; if `distribution = "NegBinom"`, the NB size is initialized from the baseline fit.
#'
#' @return An `excodeModel` whose emission is a `MultiState` formula with `nStates = states`, initialized transition matrix and `initial_mu` (and NB size when applicable).
#'
#' @seealso \code{\link{excodeFormula}}, \code{\link{excodeFamily}}, \code{\link{excodeModel}}, \code{surveillance::anscombe.residuals}, \code{surveillance::algo.farrington.assign.weights}
#'
#' @importFrom stats glm quasipoisson predict kmeans as.formula median
#' @importFrom surveillance anscombe.residuals algo.farrington.assign.weights
init_excode <- function(surv_ts, timepoint, time_units_back, distribution,
                        states, time_trend = "None", periodic_model, period_length = 52,
                        intercept = TRUE,
                        covariate_df = NULL,
                        weights_threshold = 2.58,
                        weights_threshold_baseline = 1,
                        set_baseline_state = FALSE) {
  # surv_ts
  if (!is.data.frame(surv_ts)) {
    stop("'surv_ts' must be a data.frame.")
  }
  if (nrow(surv_ts) == 0) {
    stop("'surv_ts' must not be empty.")
  }


  if (!("state" %in% names(surv_ts))) {
    surv_ts$state <- as.numeric(NA)
  }

  if (!is.numeric(surv_ts$state)) {
    surv_ts$state <- as.numeric(surv_ts$state)
  }
  # must be numeric, no coercions performed
  if (!is.numeric(surv_ts$state)) {
    stop("`surv_ts$state` must be numeric and contain only 0 or NA.")
  }
  # allow NA; all non-NA must be exactly 0
  bad_state <- !is.na(surv_ts$state) & surv_ts$state != 0
  if (any(bad_state)) {
    stop(
      "`surv_ts$state` must contain only 0 or NA (found: ",
      paste(unique(surv_ts$state[bad_state]), collapse = ", "), ")."
    )
  }

  # timepoint
  if (!is.numeric(timepoint) || length(timepoint) != 1 || timepoint <= 0) {
    stop("'timepoint' must be a positive integer index.")
  }
  if (timepoint > nrow(surv_ts)) {
    stop("'timepoint' cannot exceed the number of rows in 'surv_ts'.")
  }

  # time_units_back
  if (!is.numeric(time_units_back) || length(time_units_back) != 1 || time_units_back <= 0) {
    stop("'time_units_back' must be a positive integer.")
  }
  if (time_units_back > nrow(surv_ts)) {
    warning("'time_units_back' is larger than available rows in 'surv_ts'. Using entire dataset up to 'timepoint'.")
  }

  # distribution
  valid_distributions <- c("Poisson", "NegBinom")
  if (!distribution %in% valid_distributions) {
    stop(sprintf(
      "'distribution' must be one of: %s.",
      paste(valid_distributions, collapse = ", ")
    ))
  }

  # states
  if (!is.numeric(states) || length(states) != 1 || states < 2) {
    stop("'states' must be an integer >= 2.")
  }

  # time_trend
  valid_trends <- c("Linear", "Spline1", "Spline2", "Spline3", "Spline4", "None")
  if (!time_trend %in% valid_trends) {
    stop(sprintf(
      "'time_trend' must be one of: %s.",
      paste(valid_trends, collapse = ", ")
    ))
  }

  # periodic_model
  valid_models <- c("Mean", "FarringtonNoufaily", "Harmonic", "Custom")
  if (!periodic_model %in% valid_models) {
    stop(sprintf(
      "'periodic_model' must be one of: %s.",
      paste(valid_models, collapse = ", ")
    ))
  }

  # period_length
  if (!is.numeric(period_length) || length(period_length) != 1 || period_length <= 0) {
    stop("'period_length' must be a positive integer.")
  }

  # intercept
  if (!is.logical(intercept) || length(intercept) != 1) {
    stop("'intercept' must be TRUE or FALSE.")
  }

  # weights_threshold
  if (!is.numeric(weights_threshold) || weights_threshold < 0) {
    stop("'weights_threshold' must be >= 0.")
  }

  # weights_threshold_baseline
  if (!is.numeric(weights_threshold_baseline) || weights_threshold_baseline < 0) {
    stop("'weights_threshold_baseline' must be >= 0.")
  }


  excode_family <- excodeFamily(distribution)
  init_time_trend <- ifelse(time_trend != "None", "Linear", "None")
  excode_formula <- excodeFormula(periodic_model,
    time_trend = init_time_trend,
    intercept = intercept, data = covariate_df,
    timepoints_per_unit = period_length
  )
  excode_model <- excodeModel(
    excode_family,
    excode_formula
  )

  model_data <- extractModelData(surv_ts, excode_formula, timepoint, time_units_back)

  formula_baseline <- excode_model@emission@excode_formula@formula_bckg

  suppressWarnings(init_model <- glm(formula_baseline,
    data = model_data,
    family = quasipoisson(link = "log")
  ))
  phi <- max(summary(init_model)$dispersion, 1)
  init_anscombe_res <- surveillance::anscombe.residuals(init_model, phi)
  init_weights <- surveillance::algo.farrington.assign.weights(
    init_anscombe_res,
    weights_threshold
  )
  suppressWarnings(init_model <- glm(as.formula(formula_baseline),
    data = model_data,
    family = quasipoisson(link = "log"),
    weights = init_weights
  ))
  init_model_z <- init_model
  phi <- max(summary(init_model)$dispersion, 1)
  init_anscombe_res <- surveillance::anscombe.residuals(init_model, phi)

  has_offset <- any(names(surv_ts) == "offset")
  model_data_init <- model_data

  rtimepoints <- which(1:nrow(surv_ts) %in% model_data_init$rtime)
  surv_ts_init <- surv_ts[rtimepoints, ]

  excode_formula <- excodeFormula(periodic_model,
    time_trend = time_trend,
    intercept = intercept, data = covariate_df,
    timepoints_per_unit = period_length
  )
  excode_model <- excodeModel(
    excode_family,
    excode_formula
  )
  # model_data <- extractModelData(surv_ts, excode_formula, nrow(surv_ts), time_units_back)
  model_data <- extractModelData(surv_ts_init, excode_formula, which(timepoint == rtimepoints), time_units_back)
  params <- excode_model@emission@excode_formula@params # [-1]
  if (periodic_model != "Custom") {
    params <- params[-1]
  }
  params <- params[params != "offset(log(population))"]

  excode_formula_multistate <- excodeFormula("MultiState",
    nStates = states,
    intercept = intercept,
    time_trend = time_trend,
    offset = has_offset,
    data = model_data[, params, drop = FALSE],
    surv_ts = surv_ts_init,
    timepoints_per_unit = period_length
  )

  no_baseline_t <- which(init_anscombe_res > weights_threshold_baseline)
  baseline_t <- setdiff(1:length(init_anscombe_res), no_baseline_t)
  km <- kmeans(init_anscombe_res[no_baseline_t], states - 1)
  ordered <- order(km$centers)
  new_labels <- match(km$cluster, ordered)
  cluster <- new_labels

  initial_states <- rep(0, length(init_anscombe_res))
  for (i in 2:states) {
    initial_states[no_baseline_t[cluster == (i - 1)]] <- i - 1
  }

  model_data_init$state <- factor(initial_states)
  transMat <- table(
    model_data_init$state[1:(nrow(model_data_init) - 1)],
    model_data_init$state[2:(nrow(model_data_init))]
  ) + max(c(1, nrow(model_data_init) * 0.05))
  transMat <- t(apply(transMat, 1, function(x) x / sum(x)))
  suppressWarnings(init_model <- glm(paste0(formula_baseline, " + state"),
    data = model_data_init,
    family = quasipoisson(link = "log")
  ))
  phi <- max(summary(init_model)$dispersion, 1)

  initial_mu <- matrix(NA,
    nrow = nrow(model_data_init),
    ncol = states
  )
  for (i in 1:states) {
    newdata <- model_data_init
    newdata$state <- factor(i - 1, levels = 0:(states - 1))
    initial_mu[, i] <- predict(init_model, newdata = newdata, type = "response")
  }


  theta <- median(predict(init_model, type = "response") / (phi - 1))
  excode_family <- excodeFamily(distribution, nb_size = theta)
  excode_multistate <- excodeModel(excode_family,
    excode_formula_multistate,
    initial_mu = initial_mu,
    transMat = transMat
  )

  if (set_baseline_state) {
    excode_multistate <- compute_baseline_state(
      excode_multistate,
      weights_threshold_baseline
    )
  }

  excode_multistate
}


#' Update EXCODE model by setting baseline states in surv_ts
#'
#' @description
#' Fits a quasi-Poisson GLM using the **background formula** stored in
#' \code{excode_multistate@emission@excode_formula@formula_bckg} on the data found in
#' \code{excode_multistate@emission@excode_formula@surv_ts} (optionally augmented with
#' covariates from \code{@data} if present), computes Anscombe residuals, and **writes**
#' a \code{state} column back into \code{@surv_ts}: entries are \code{0} when residuals
#' are below \code{weights_threshold_baseline} (or the observed count is \code{0}), and
#' \code{NA} otherwise. Returns the updated \code{excodeModel}.
#'
#' @param excode_multistate An \code{excodeModel} returned by \code{init_excode()}
#'   with a \code{"MultiState"} emission.
#' @param weights_threshold_baseline Numeric residual cutoff; residuals strictly below
#'   this threshold are labeled baseline (\code{0}). Default \code{1}.
#'
#' @return The **updated** \code{excodeModel} with
#'   \code{excode_multistate@emission@excode_formula@surv_ts$state} set.
#'
#' @importFrom stats glm quasipoisson
#' @importFrom surveillance anscombe.residuals
#' @export
compute_baseline_state <- function(excode_multistate,
                                   weights_threshold_baseline = 1) {
  # --- validation ---
  if (missing(excode_multistate) || is.null(excode_multistate)) {
    stop("'excode_multistate' must be provided and must be an excodeModel.")
  }
  if (is.null(excode_multistate@emission) ||
    is.null(excode_multistate@emission@excode_formula)) {
    stop("'excode_multistate' lacks a valid emission/excode_formula slot.")
  }
  if (!is.numeric(weights_threshold_baseline) ||
    length(weights_threshold_baseline) != 1 ||
    !is.finite(weights_threshold_baseline) ||
    is.na(weights_threshold_baseline)) {
    stop("'weights_threshold_baseline' must be a finite number.")
  }

  xf <- excode_multistate@emission@excode_formula

  # Pull surv_ts (required) and covariates data (optional)
  surv_ts <- tryCatch(xf@surv_ts, error = function(e) NULL)
  if (is.null(surv_ts) || !is.data.frame(surv_ts) || nrow(surv_ts) == 0) {
    stop("No valid data in 'excode_multistate@emission@excode_formula@surv_ts'.")
  }
  data_cov <- tryCatch(xf@data, error = function(e) NULL)

  # Combine covariates from @data if present (append missing columns by row order)
  model_df <- surv_ts
  if (!is.null(data_cov)) {
    if (!is.data.frame(data_cov)) {
      warning("@data inside excode_formula is not a data.frame; ignoring it.")
    } else if (nrow(data_cov) != nrow(model_df)) {
      stop("Row count mismatch between '@surv_ts' and '@data' in the model.")
    } else {
      extra_cols <- setdiff(names(data_cov), names(model_df))
      if (length(extra_cols)) {
        model_df <- cbind(model_df, data_cov[, extra_cols, drop = FALSE])
      }
    }
  }

  # Background formula for GLM
  if (is.null(xf@formula_bckg)) {
    stop("No background formula '@formula_bckg' in the model.")
  }
  formula_bckg <- xf@formula_bckg

  # Identify response in @surv_ts (you indicated it's 'observed')
  response_col <- "observed"
  if (!response_col %in% names(model_df)) {
    stop("Expected a numeric 'observed' column in '@surv_ts' for the response.")
  }
  if (!is.numeric(model_df[[response_col]])) {
    stop("'observed' in '@surv_ts' must be numeric.")
  }
  if (any(model_df[[response_col]] < 0, na.rm = TRUE)) {
    stop("'observed' in '@surv_ts' contains negative values.")
  }

  # Ensure a 'population' column exists if the formula references offset(log(population))
  if (!"population" %in% names(model_df)) {
    model_df$population <- 1
  }

  # Some pipelines store 'response' in the formula; map it to 'observed' if needed
  if (is.character(formula_bckg)) {
    formula_bckg <- gsub("\\bresponse\\b", "observed", formula_bckg)
  }

  # Fit GLM using the model's background formula on the constructed data
  suppressWarnings(
    curr_model <- stats::glm(formula_bckg,
      data = model_df,
      family = stats::quasipoisson(link = "log")
    )
  )

  # Dispersion and residuals
  phi <- max(summary(curr_model)$dispersion, 1)
  curr_anscombe_res <- surveillance::anscombe.residuals(curr_model, phi)
  if (length(curr_anscombe_res) != nrow(model_df)) {
    stop("Internal error: residual vector length does not match data rows.")
  }

  # Build baseline state labels: 0 if residual < threshold OR observed == 0; else NA
  state_vec <- rep(NA_integer_, nrow(model_df))
  obs <- model_df[[response_col]]
  set_zero_idx <- which((curr_anscombe_res < weights_threshold_baseline | obs == 0) &
    !is.na(curr_anscombe_res) & !is.na(obs))
  if (length(set_zero_idx)) state_vec[set_zero_idx] <- 0L

  # Helpful edge-case warnings
  if (all(is.na(state_vec))) {
    warning("All states are NA; no observations passed the baseline threshold rule.")
  } else if (all(!is.na(state_vec) & state_vec == 0L)) {
    warning("All states are 0; threshold may be too high or data may be zero-inflated.")
  }

  # If the MultiState model has nStates, set factor levels accordingly; otherwise leave integer
  n_states <- tryCatch(xf@nStates, error = function(e) NULL)
  surv_ts$state <- as.numeric(state_vec)
  surv_ts$state[nrow(surv_ts)] <- NA

  # Write updated surv_ts back to the model and return it
  excode_multistate@emission@excode_formula@surv_ts <- surv_ts
  excode_multistate
}
