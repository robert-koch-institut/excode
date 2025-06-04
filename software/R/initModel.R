# Initialization procedure for semi- and unsupervised learning
initGLM <- function(hmm,
                    modelData,
                    learning = "unsupervised",
                    weightsThreshold = 2.58,
                    weightsThresholdBckg = 1,
                    setBckgState = TRUE,
                    limitCaseWeeks = c(-Inf, 1)) {
  rownames(modelData) <- NULL
  suppressWarnings(curr_model <- glm(hmm@emission@excode_formula@formula_bckg,
    data = modelData[modelData$init, ],
    family = quasipoisson(link = "log")
  ))

  phi <- max(summary(curr_model)$dispersion, 1)
  K <- min(c(10, phi + 1))

  curr_anscombe_res_temp <- surveillance::anscombe.residuals(curr_model, max(summary(curr_model)$dispersion, 1))
  curr_anscombe_res <- rep(Inf, nrow(modelData))
  curr_anscombe_res[which(modelData$init)] <- curr_anscombe_res_temp # & bckg_state

  # Use anscombe residuals to set bckg states
  modelData$bckg_state <- modelData$true_state

  modelData$bckg_state <- NA
  if (setBckgState) {
    selInd <- which(!modelData$curr_week)
    setZero <- which(curr_anscombe_res[selInd] < weightsThresholdBckg | modelData$response[selInd] == 0)
    modelData$bckg_state[setZero] <- 0
  }
  # Set states which should not be used for training to NA
  modelData$true_state[!modelData$state_training] <- NA
  modelData$bckg_state[!modelData$state_training] <- NA

  if (learning == "unsupervised" | (learning == "semisupervised" & all(is.na(modelData$true_state)))) {
    modelData$known_state <- modelData$bckg_state
  } else {
    modelData$known_state <- modelData$true_state
  }
  # write.table(modelData, file="data_excode.csv", sep=";", quote=FALSE, row.names=FALSE)

  dataGLM1 <- modelData
  dataGLM2 <- modelData
  dataGLM1$state <- 0
  dataGLM2$state <- 1
  dataGLM_expanded <- rbind(dataGLM1, dataGLM2)

  expected_mu_bckg <- predict(curr_model, newdata = modelData, type = "response")
  expected_mu <- c(expected_mu_bckg, expected_mu_bckg * K)
  start_pi <- sum(as.numeric(dataGLM_expanded$response == 0)) / (nrow(dataGLM_expanded))
  emission_params <- data.frame(
    mu = expected_mu, nb_size = expected_mu / (phi - 1),
    pi = start_pi
  )
  take_cols <- names(dataGLM_expanded)

  # Call generic function for Negative Binomial for initialization
  emissionProb <- calcEmissionProb(new("NegBinom"), cbind(dataGLM_expanded, emission_params))

  if (all(emissionProb[, 2] == 0)) {
    K <- max(c(1.2, max(modelData$response / expected_mu_bckg)))
    expected_mu_bckg <- predict(curr_model, newdata = modelData, type = "response")
    expected_mu <- c(expected_mu_bckg, expected_mu_bckg * K)
    start_pi <- sum(as.numeric(dataGLM_expanded$response == 0)) / (nrow(dataGLM_expanded))
    emission_params <- data.frame(
      mu = expected_mu, nb_size = expected_mu / (phi - 1),
      pi = start_pi
    )
    take_cols <- names(dataGLM_expanded)
    emissionProb <- calcEmissionProb(new("NegBinom"), cbind(dataGLM_expanded, emission_params))
  }

  dataGLM_expanded <- cbind(dataGLM_expanded, emission_params)

  dataGLM_expanded <- dataGLM_expanded[, c(take_cols, getParamNames(hmm@emission@distribution))]

  both_zero <- which(apply(emissionProb, 1, sum) < .Machine$double.eps)
  if (length(both_zero) > 0) {
    for (k in both_zero) {
      assign_state <- ifelse(modelData$response[k] > expected_mu[k], 1, 0)
      emissionProb[k, ] <- c(1 - assign_state, assign_state)
    }
  }


  # gamma_xsi = getGammaXsi(hmm, dataGLM_expanded, emissionProb)
  gamma_xsi <- forwardBackward(
    hmm, dataGLM_expanded,
    emissionProb
  )

  # write.table(gamma_xsi$gamma, file="posterior_excode.csv", sep=";", quote=FALSE, row.names=FALSE)

  omega <- as.vector(gamma_xsi$gamma)
  model_updated <- updateEmission(hmm@emission, dataGLM_expanded, omega)
  model_updated
}

init_glm_mutlistate <- function(hmm,
                                modelData,
                                learning = "unsupervised",
                                setBckgState = TRUE) {
  nStates <- hmm@nStates
  ts_len <- nrow(modelData)
  dataGLM_expanded <- modelData[rep(1:nrow(modelData), nStates), ]
  subset_mu <- modelData$rtime + (as.numeric(as.factor(modelData$id)) - 1) * max(modelData$rtime)
  dataGLM_expanded$mu <- as.vector(hmm@emission@mu[subset_mu, ])
  dataGLM_expanded$state <- sort(rep(0:(nStates - 1), nrow(modelData)))

  # start_state <- ifelse(hmm@emission@excode_formula@intercept, 0, 1)
  for (curr_state in 1:(nStates - 1)) {
    dataGLM_expanded[[paste0("state", curr_state)]] <- as.numeric(curr_state == dataGLM_expanded$state)
  }
  # dataGLM_expanded %>% View()
  # stop("nonono")
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
    if (length(nb_size) == length(unique(modelData$id))) {
      nb_size <- unlist(lapply(nb_size, function(x) rep(x, length(unique(modelData$rtime)))))
      nb_size <- rep(nb_size, nStates)
    } else {
      nb_size <- nb_size[subset_mu]
    }
    dataGLM_expanded$nb_size <- nb_size
  }

  if (learning == "unsupervised" | (learning == "semisupervised" &
    all(is.na(dataGLM_expanded$true_state)))) {
    dataGLM_expanded$known_state <- dataGLM_expanded$true_state
  } else {
    dataGLM_expanded$known_state <- dataGLM_expanded$true_state
  }


  emissionProb <- calcEmissionProb(
    hmm@emission@distribution,
    dataGLM_expanded
  )
  # View(cbind(emissionProb, matrix(dataGLM_expanded$mu, ncol=3), matrix(dataGLM_expanded$response, ncol=3)))
  # stop("a")
  # state_seq <- apply(emissionProb,1,function(x) which(x==max(x))[1])
  # transMat <- table(state_seq[1:(length(state_seq)-1)], state_seq[2:length(state_seq)])+1
  # transMat <- t(apply(transMat,1,function(x) x/sum(x)))
  # print(apply(transMat,1,sum))
  # hmm@transitions <- transMat
  gamma_xsi <- forwardBackward(
    hmm, dataGLM_expanded,
    emissionProb
  )
  # print(round(gamma_xsi$gamma,4))
  omega <- as.vector(gamma_xsi$gamma)
  model_updated <- updateEmission(hmm@emission, dataGLM_expanded, omega)
  model_updated
}
