setMethod("updateEmission",
  signature = c("EmissionGLMNegBinom", "data.frame", "numeric"),
  function(emission, dat, omega) {
    formula <- emission@excode_formula@formula
    non_zero_weights <- which(omega != 0)
    omega_temp <- omega[non_zero_weights]
    dat_orig <- dat
    dat <- dat[non_zero_weights, ]

    suppressWarnings(fit <- glm.nb(as.formula(formula), dat,
      weights = omega_temp
    ))

    mu <- predict(fit, newdata = dat_orig, type = "response")
    emission@excode_formula@coefficients <- fit$coefficients

    dat_orig$mu <- mu
    dat_orig$nb_size <- fit$theta
    class(fit) <- c("glm", "lm")
    emission@glm <- fit
    list(model = dat_orig, emission = emission)
  }
)
