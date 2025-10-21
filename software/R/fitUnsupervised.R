fitUnsupervised <- function(hmm, modelData, transMat_init, maxIter, verbose, time_units_back) {
  model_init <- NULL
  model <- NULL

  model_init <- init_glm_mutlistate(hmm, modelData,
    setBckgState = hmm@setBckgState
  )

  model <- model_init$model
  hmm@emission <- model_init$emission
  base_weight <- matrix(
    min(c(
      nrow(modelData),
      time_units_back * hmm@emission@excode_formula@timepoints_per_unit
    )),
    nrow = hmm@nStates, ncol = hmm@nStates
  )

  # Calculate stationary state distribution
  n <- ncol(transMat_init)
  A <- t(transMat_init - diag(n))
  A <- rbind(A, rep(1, n))
  b <- c(rep(0, n), 1)
  stationary_state <- qr.solve(A, b)
  # Calculate prior weights
  prior_weights <- round(transMat_init * base_weight * stationary_state + 1)
  hmm@prior_weights <- prior_weights
  initProb <- hmm@initial_prob
  hmm@transitions <- matrix(1 / hmm@nStates, nrow = hmm@nStates, ncol = hmm@nStates)
  emission_prob <- calcEmissionProb(hmm@emission@distribution, model)
  hmm_expectation <- forwardBackward(hmm, model, emission_prob)

  old_loglik <- -Inf

  if (hmm@transitions_prior) {
    log_dir_ll <- 0
    for (i in 1:nrow(hmm@transitions)) {
      log_dir_ll <- log_dir_ll + log_dirichlet(
        hmm@transitions[i, ],
        hmm@prior_weights[i, ]
      )
    }
    hmm@loglik_transitions <- log_dir_ll # log_dir1 + log_dir2
  }

  new_loglik <- hmm_expectation$LogLik + hmm@loglik_transitions
  curr_diff <- Inf

  if (verbose) {
    cat(" Iter 0: LogLik=", new_loglik, " (diff=", curr_diff, ")\n", sep = "")
  }

  niter <- 1

  # EM iterations
  while (curr_diff > 1e-6 & niter <= maxIter) {
    old_loglik <- new_loglik
    omega <- as.vector(hmm_expectation$gamma)
    model_updated <- updateEmission(hmm@emission, model, omega)
    hmm@emission <- model_updated$emission
    model <- model_updated$model
    hmm <- updateTransInitProb(hmm, hmm_expectation$gamma, hmm_expectation$xsi, model)

    emission_prob <- calcEmissionProb(hmm@emission@distribution, model)
    hmm_expectation <- forwardBackward(hmm, model, emission_prob)
    new_loglik <- hmm_expectation$LogLik + hmm@loglik_transitions

    curr_diff <- (new_loglik - old_loglik) / abs(old_loglik)
    if (verbose) {
      cat(" Iter ", niter, ": LogLik=", new_loglik, " (diff=", curr_diff, ")\n", sep = "")
    }
    niter <- niter + 1
  }

  if (niter > maxIter) {
    warning("Maximum number of iterations exceeded in EM.")
  }

  if (curr_diff < 0) {
    stop("Decrease in Log-likelihood during EM.\n")
  }

  list(
    hmm = hmm, hmm_expectation = hmm_expectation, model = model,
    loglik = new_loglik, niter = niter
  )
}
