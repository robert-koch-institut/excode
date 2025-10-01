setMethod("updateEmission",
  signature = c("EmissionGLMPoisson", "data.frame", "numeric"),
  function(emission, dat, omega) {
    formula <- emission@excode_formula@formula
    mf <- model.frame(as.formula(formula), data = dat)

    # Extract response, model matrix and offset
    model_offset_all <- model.offset(mf)
    y <- model.response(mf, type = "numeric")
    modelData <- model.matrix(as.formula(formula), dat)
    X <- modelData
    # Remove data with weight==0 for faster and more stable fitting
    non_zero_weights <- which(omega != 0)
    modelData <- modelData[non_zero_weights, ]
    model_offset <- model_offset_all[non_zero_weights]
    y <- y[non_zero_weights]
    omega_temp <- omega[non_zero_weights]
    dat_orig <- dat
    dat <- dat[non_zero_weights, ]

    etastart <- NULL
    if (FALSE & !any(is.na(emission@excode_formula@coefficients))) {
      etastart <- modelData %*% emission@excode_formula@coefficients + model_offset
    }
    # (modelData %>% View())
    # print(all(modelData==0))
    # files <- list.files("./test_emission/")
    suppressWarnings(fit <- glm.fit(
      x = modelData, y = y,
      etastart = etastart,
      family = poisson(link = "log"),
      weights = omega_temp,
      offset = model_offset
    ))
    emission@excode_formula@coefficients <- fit$coefficients
    # Calclulate updated mean and nb_size values for all time series
    mu <- exp(X %*% fit$coefficients + model_offset_all)
    # print(fit$coefficients)
    # save(fit, mu, modelData, y, etastart,omega_temp, non_zero_weights, omega, dat,
    #     file=paste0("./test_emission/test_", length(files)+1, ".rda"))

    dat_orig$mu <- as.vector(mu)
    # View(dat_orig)
    # stop("ss")
    list(model = dat_orig, emission = emission)
  }
)
