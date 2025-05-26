calcPosteriorBound <- function(excode_model, min_posterior, maxiter = 1000) {
  result_list <- list()
  ts_id <- unique(excode_model@id)

  for (timepoint_fit in unique(excode_model@alpha[, 4])) {
    for (id in ts_id) {
      alpha_index <- which(excode_model@alpha[, 4] == timepoint_fit &
        rownames(excode_model@alpha) == id)
      alpha_init <- excode_model@alpha[alpha_index, ]
      alpha_timepoint <- alpha_init[, 3]
      alpha_init <- alpha_init[, 1:2]
      transMat <- excode_model@transitions

      index <- which(excode_model@id == id &
        excode_model@timepoint_fit == timepoint_fit &
        excode_model@timepoint %in% alpha_timepoint)

      mu_index <- which(excode_model@timepoint == timepoint_fit &
        excode_model@timepoint_fit == timepoint_fit &
        excode_model@id == id)

      mu0 <- excode_model@emission@mu0[mu_index]
      mu1 <- excode_model@emission@mu1[mu_index]

      alpha0 <- c()
      alpha1 <- c()

      min_posterior_temp <- -Inf
      y <- excode_model@observed[mu_index]
      if (alpha_init[2, 2] >= min_posterior) {
        min_posterior_temp <- 1 - min_posterior
        alpha_cutoff <- alpha_init[2, 1]
      } else {
        min_posterior_temp <- min_posterior
        alpha_cutoff <- alpha_init[2, 2]
      }

      niter <- 0
      while (alpha_cutoff[length(alpha_cutoff)] < min_posterior_temp & niter <= maxiter) {
        alpha <- alpha_init[rownames(alpha_init) == id, ]

        curr_data <- create_emission_prob_input(
          excode_model@emission@distribution,
          y, mu0, mu1, index
        )

        emissionProb <- calcEmissionProb(
          excode_model@emission@distribution,
          curr_data
        )

        t <- 2
        resc <- 0
        for (i in 1:nrow(transMat)) {
          alpha[t, i] <- 0
          for (j in 1:ncol(transMat)) {
            alpha[t, i] <- alpha[t, i] + alpha[t - 1, j] * transMat[j, i]
          }
          alpha[t, i] <- alpha[t, i] * emissionProb[t - 1, i]
          resc <- resc + alpha[t, i]
        }
        # Rescale alpha[t][i]
        resc <- 1 / resc

        for (i in 1:nrow(transMat)) {
          alpha[t, i] <- resc * alpha[t, i]
        }
        alpha0[length(alpha0) + 1] <- alpha[t, 1]
        alpha1[length(alpha0) + 1] <- alpha[t, 2]

        if (alpha_init[2, 2] >= min_posterior) {
          y <- y - 1
          alpha_cutoff[length(alpha_cutoff) + 1] <- alpha[t, 1]
        } else {
          y <- y + 1
          alpha_cutoff[length(alpha_cutoff) + 1] <- alpha[t, 2]
        }
        niter <- niter + 1
      }

      posterior_ub <- NA
      if (niter == maxiter) {
        posterior_ub <- NA
      }

      if (!all(is.na(alpha1))) {
        if (alpha_init[2, 2] >= min_posterior) {
          posterior_ub <- y + 2
        } else {
          posterior_ub <- y - 1
        }
      }


      result_list[[length(result_list) + 1]] <-
        data.frame(
          posterior_ub = posterior_ub,
          timepoint_fit = timepoint_fit,
          id = id
        )
    }
  }
  result <- do.call("rbind", result_list)
  rownames(result) <- NULL
  result
}




setGeneric("calcPvalueBound", function(distribution, excode_model, alpha) standardGeneric("calcPvalueBound"))

setMethod("calcPvalueBound",
  signature = c("NegBinom", "excodeModel", "numeric"),
  function(distribution, excode_model, alpha) {
    index <- which(excode_model@timepoint_fit == excode_model@timepoint)
    ub_pval <- qnbinom(alpha,
      size = excode_model@emission@distribution@nb_size[index],
      mu = excode_model@emission@mu0[index],
      lower.tail = FALSE
    )
    data.frame(
      pval_ub = ub_pval,
      timepoint_fit = excode_model@timepoint_fit[index],
      id = excode_model@id[index]
    )
  }
)

setMethod("calcPvalueBound",
  signature = c("Poisson", "excodeModel", "numeric"),
  function(distribution, excode_model, alpha) {
    index <- which(excode_model@timepoint_fit == excode_model@timepoint)

    ub_pval <- qpois(alpha, lambda = excode_model@emission@mu0[index], lower.tail = FALSE)
    data.frame(
      pval_ub = ub_pval,
      timepoint_fit = excode_model@timepoint_fit[index],
      id = excode_model@id[index]
    )
  }
)
