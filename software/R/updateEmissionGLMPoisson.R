setMethod("updateEmission",
  signature = c("EmissionGLMPoisson", "data.frame", "numeric"),
  function(emission, dat, omega) {
    formula <- emission@excode_formula@formula
    non_zero_weights <- which(omega != 0)
    omega_temp <- omega[non_zero_weights]
    dat_orig <- dat
    dat <- dat[non_zero_weights, ]

    fit <- glm(as.formula(formula), dat,
      family = quasipoisson(),
      weights = omega_temp
    )

    mu <- predict(fit, newdata = dat_orig, type = "response")
    emission@excode_formula@coefficients <- fit$coefficients
    emission@glm <- fit

    dat_orig$mu <- mu
    list(model = dat_orig, emission = emission)
  }
)
