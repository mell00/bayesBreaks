
baar <- function(k = NULL, time, data, iterations, burn_in = 50, make_murder_p = 0.5,
  percent = 0.02, lambda = 1, jump_p = 0.25, ar = 1, progress = TRUE, fit_storage = FALSE) {

  ar <- floor(ar)

  if (length(time) != length(data)) {
    stop("Data and time vectors must be of equal length.")
  } else if (length(data) < (6 * ar)) {
    stop("Data insufficient for order of AR model. Try a lower order.")
  } else if (make_murder_p >= 1) {
    stop("Make/murder proportion must be less than 1.")
  } else if (percent >= 0.5) {
    stop("Percent for jiggle neighborhood must be less than 0.5.")
  } else if (jump_p > 1) {
    stop("Jump proportion must be less than or equal to 1.")
  }

  full_data <- cbind(seq_along(as.numeric(time)), as.numeric(data))  # combining time and data inputs
  n <- length(full_data[, 1])  # number of observations
  k_ends <- c(min(full_data[, 1]), stats::na.omit(k), n)  # adding end points to k

  extract_residuals <- function(model) {
    res <- NULL
    if (!is.null(model$res)) {
      res <- model$res
    } else if (!is.null(model$resid)) {
      res <- model$resid
    } else {
      res <- tryCatch(stats::residuals(model), error = function(e) NULL)
    }
    if (is.null(res)) {
      stop("Unable to extract residuals from fitted AR model.")
    }
    return(res)
  }

  fit_segment <- function(values) {
    suppressWarnings(FitAR(values, p = ar, demean = TRUE))
  }

  # Compute the sum of log likelihoods across all segments.
  fitMetrics <- function(k_ends, full_data) {

    sum_loglik <- 0

    if (length(k_ends) < 3) {
      model <- fit_segment(full_data[, 2])
      see <- sum(stats::na.omit(extract_residuals(model))^2)
      s2 <- see/n
      sum_loglik <- (-1 * n/2) * (log(2 * pi) + log(s2) + 1)
    } else {
      for (i in 2:length(k_ends)) {
        start_idx <- if (i == 2)
          k_ends[i - 1] else k_ends[i - 1] + 1
        end_idx <- k_ends[i]
        y_values <- full_data[start_idx:end_idx, 2]
        model <- fit_segment(y_values)
        sub_n <- length(y_values)
        see <- sum(stats::na.omit(extract_residuals(model))^2)
        s2 <- see/sub_n
        sub_loglik <- (-1 * sub_n/2) * (log(2 * pi) + log(s2) + 1)
        sum_loglik <- sum_loglik + sub_loglik
      }
    }
    sum_loglik
  }

  # Identify all available positions for placing new breakpoints.
  freeObservations <- function(k_ends, ar) {
    n <- max(k_ends)
    full_set <- seq_len(n)
    constraint <- if (ar == 1) {
      2
    } else {
      (2 * ar - 1)
    }

    neighborhoods <- lapply(k_ends, function(idx) {
      seq.int(max(1, idx - constraint), min(n, idx + constraint))
    })
    exclude_set <- unique(unlist(neighborhoods, use.names = FALSE))
    diff_set <- setdiff(full_set, exclude_set)
    return(sort(diff_set))
  }

  # Randomly add a breakpoint.
  barMake <- function(k_ends) {

    diff_set <- freeObservations(k_ends, ar)
    if (length(diff_set) < 2) {
      rand_spot <- diff_set
    } else {
      rand_spot <- sample(diff_set, 1)
    }
    k_ends_final <- sort(c(k_ends, rand_spot))
    return(k_ends_final)

  }


  # Randomly remove a breakpoint.
  barMurder <- function(k_ends) {

    k <- k_ends[c(-1, -length(k_ends))]
    if (length(k) == 1) {
      random_num <- 1
    } else {
      random_num <- sample(1:length(k), 1)
    }
    k_ends_final <- k_ends[-(random_num + 1)]
    return(k_ends_final)

  }

  # Relocate a breakpoint without respecting neighbourhood constraints.
  barMove <- function(k_ends) {

    k_ends_less <- barMurder(k_ends)
    k_ends_final <- barMake(k_ends_less)
    return(k_ends_final)

  }

  # Move a breakpoint within a local neighbourhood.
  barJiggle <- function(percent, k_ends) {

    k <- k_ends[c(-1, -length(k_ends))]
    random_num <- sample(1:length(k), 1)
    random_bkpt <- k_ends[random_num + 1]
    full_data <- c(1:max(k_ends))
    wiggliness <- max(k_ends) * percent
    if (floor(random_bkpt - wiggliness) > 0) {
      left_lim <- floor(random_bkpt - wiggliness)
    } else {
      left_lim <- 1
    }
    if (ceiling(random_bkpt + wiggliness) < max(k_ends)) {
      right_lim <- ceiling(random_bkpt + wiggliness)
    } else {
      right_lim <- max(k_ends)
    }
    prelim_neighborhood <- full_data[left_lim:right_lim]
    left_neighbor <- k_ends[random_num]
    right_neighbor <- k_ends[random_num + 2]
    if (ar == 1) {
      constraint <- (3 - 1)
    } else {
      constraint <- (2 * ar - 1)
    }
    lr_limit <- left_neighbor + constraint
    rl_limit <- right_neighbor - constraint
    exclusions <- sort(c(full_data[1:lr_limit], full_data[rl_limit:max(k_ends)]))
    final_neighborhood <- setdiff(prelim_neighborhood, exclusions)

    if (length(final_neighborhood) > 1) {
      new_location <- sample(final_neighborhood, 1)
      k_ends_less <- k_ends[-(random_num + 1)]
      final_k_ends <- sort(c(k_ends_less, new_location))
      return(final_k_ends)
    } else if (length(final_neighborhood) == 1) {
      new_location <- final_neighborhood
      k_ends_less <- k_ends[-(random_num + 1)]
      final_k_ends <- sort(c(k_ends_less, new_location))
      return(final_k_ends)
    } else {
      return("jiggle failure")
    }

  }

  # Propose a new breakpoint set and report proposal probabilities and step
  # type.
  newEnds <- function(k_ends, make_k, murder_k) {

    u_step <- runif(1)  # random number from 0 to 1 taken from a uniform distribution for selecting step

    if (ar == 1) {
      constraint <- 5
    } else {
      constraint <- ar * 4
    }

    if (max(diff(k_ends)) >= constraint & length(k_ends) < 3 | max(diff(k_ends)) >=
      constraint & u_step <= make_k) {
      type <- "add"
      a.count <<- a.count + 1
      k_ends_new <- barMake(k_ends)  # make

      # setting up qs for ratio
      q1 <- murder_k/(length(k_ends_new) - 2)
      n_free <- length(freeObservations(k_ends, ar))
      q2 <- make_k/n_free

    } else if (u_step > make_k & u_step <= (make_k + murder_k)) {
      type <- "sub"
      s.count <<- s.count + 1
      k_ends_new <- barMurder(k_ends)  # murder

      # setting up qs for ratio
      n_free <- length(freeObservations(k_ends_new, ar))
      q1 <- make_k/n_free
      q2 <- murder_k/(length(k_ends) - 2)

    } else {
      move_u <- runif(1)
      if (move_u < jump_p) {
        type <- "move"
        m.count <<- m.count + 1
        k_ends_new <- barMove(k_ends)  # jump

        # fake qs because they cancel
        q1 <- 1
        q2 <- 1
      } else {
        type <- "jiggle"
        j.count <<- j.count + 1
        k_ends_new <- barJiggle(percent, k_ends)  # jiggle

        if (k_ends_new[[1]] == "jiggle failure") {
          k_ends_new <- k_ends
        }

        # fake qs because they cancel
        q1 <- 1
        q2 <- 1
      }
    }

    return(list(k_ends_new, q1, q2, type))
  }

  # setting up counters for burn-in Metropolis-Hastings
  type <- "0"
  a.count <- 0
  s.count <- 0
  m.count <- 0
  j.count <- 0
  add.accept.count <- 0
  sub.accept.count <- 0
  move.accept.count <- 0
  jiggle.accept.count <- 0

  # getting constants for qs for burn-in Metropolis-Hastings
  starting_bkpts <- length(k_ends) - 1  # most probable number of breakpoints based on starting info
  starting_nfree <- length(freeObservations(k_ends, ar))  # most probable n_free based on starting info
  starting_ttl <- starting_bkpts + starting_nfree  # total to get percentages
  make_k <- make_murder_p * (starting_nfree/starting_ttl)  # proportion for make
  murder_k <- make_murder_p * (starting_bkpts/starting_ttl)  # proportion for murder

  if (progress == TRUE) {
    writeLines("\nBeginning burn period.")
    burn_progress <- utils::txtProgressBar(min = 0, max = burn_in, style = 3)
  }

  # Burn Metropolis Hasting
  for (i in 1:burn_in) {

    old_loglik <- fitMetrics(k_ends, full_data)  # gets log likelihood for existing breakpoints

    k_and_q <- newEnds(k_ends, make_k, murder_k)
    k_ends_new <- k_and_q[[1]]
    q1 <- k_and_q[[2]]
    q2 <- k_and_q[[3]]

    new_loglik <- fitMetrics(k_ends_new, full_data)

    # birth and death ratios


    old_new_product_birth <- function(k) {
      prod(n - (3 * ar - (q1 + 1) * (2 * ar + k)))
    }
    product_death <- function(k) {
      prod(n - 3 * ar - q1 * 2 * ar + k)
    }
    old_new_product_death <- function(k) {
      prod(n - 3 * ar - (q1 - 1) * (2 * ar) + k)
    }


    birth_old_to_new_ratio <- (factorial(starting_bkpts + 1) * q2 * (dpois(length(k_ends_new) -
      2, lambda)))/(old_new_product_birth(starting_bkpts))
    birth_new_to_old_ratio <- (factorial(starting_bkpts) * q1 * (dpois(length(k_ends) -
      2, lambda)))/(product_death(starting_bkpts))
    death_new_to_old_ratio <- (factorial(starting_bkpts) * q2 * (dpois(length(k_ends) -
      2, lambda)))/(product_death(starting_bkpts))
    death_old_to_new_ratio <- (factorial(starting_bkpts - 1) * make_k * (dpois(length(k_ends_new) -
      2, lambda)))/(old_new_product_death(starting_bkpts))
    # end of birth and death ratios

    delta_bic <- (-2 * new_loglik + log(n) * (length(k_ends_new) - 1) * (3 +
      ar)) - (-2 * old_loglik + log(n) * (length(k_ends) - 1) * (3 + ar))
    ratio <- (-1 * delta_bic/2) + log(birth_old_to_new_ratio) - log(death_new_to_old_ratio)
    u_ratio <- log(runif(1))  # random number from 0 to 1 taken from a uniform distribution and then log transformed

    if (abs(delta_bic) == Inf) {
      # safe guard against random models creating infinite ratios
      k_ends <- k_ends  # old
    } else if (ratio >= u_ratio) {
      k_ends <- k_ends_new  # new
    } else {
      k_ends <- k_ends  # old
    }

    if (progress == TRUE) {
      utils::setTxtProgressBar(burn_progress, i)
    }

  }

  # initializing matrices/storage objects for final Metropolis-Hasting
  all_k_best <- data.frame(matrix(ncol = (length(k_ends) - 2), nrow = 0))
  if (fit_storage == TRUE) {
    beta_draws <- list()
    sigma_draws <- list()
    all_fits <- matrix(numeric(0), nrow = 0, ncol = n)
    all_MSE <- numeric(0)
  }
  all_BIC <- data.frame()
  accept_count <- 0

  # setting up counters for final Metropolis-Hasting
  type <- "0"
  a.count <- 0
  s.count <- 0
  m.count <- 0
  j.count <- 0
  add.accept.count <- 0
  sub.accept.count <- 0
  move.accept.count <- 0
  jiggle.accept.count <- 0

  # setting up priors for beta draws (define what b_0 and B_0 are)
  if (fit_storage == TRUE) {
    alt_arima <- function(full_data, ar) {
      tryCatch(arima(full_data[, 2], method = "ML", order = c(ar, 0, 0)), error = function(e) arima(full_data[,
        2], method = "CSS", order = c(ar, 0, 0)))
    }
    model <- suppressWarnings(alt_arima(full_data, ar))
    informationless <- matrix(0, ncol = (ar + 1), nrow = (ar + 1))
    diag(informationless) = rep(1000, (ar + 1))
    alt_solve <- function(model_coef) {
      tryCatch(solve(model_coef), error = function(e) informationless)
    }
    fisher <- suppressWarnings(alt_solve(model$var.coef))  # amount of data contained in 1 data point
    smiley <- n * fisher  # empirical Bayes (using data to set priors) # as n goes to inf, variance becomes unbiased

    coef_list <- model$coef[[length(model$coef)]]

    for (a in 1:(length(model$coef) - 1)) {

      coef_list <- c(coef_list, model$coef[[a]], recursive = T)  # pulls each coefficient from model

    }

    b_0 <- matrix(coef_list, (ar + 1), 1)  # matrix of beta means for posterior draw
    B_0 <- smiley  # variance-covariance matrix for posterior draw

    # beta and sigma draw
    post_beta_list <- data.frame(Empty = rep(NA, (ar + 1)))
    post_sigma_list <- data.frame(Empty = NA)
    prior_precision <- tryCatch(solve(B_0), error = function(e) MASS::ginv(B_0))
  }

  # getting constants for qs for final Metropolis-Hasting
  starting_bkpts_2 <- length(k_ends) - 1  # most probable number of breakpoints based on starting info
  starting_nfree <- length(freeObservations(k_ends, ar))
  starting_ttl <- starting_bkpts_2 + starting_nfree  # total to get percentages
  make_k <- make_murder_p * (starting_nfree/starting_ttl)  # proportion for make
  murder_k <- make_murder_p * (starting_bkpts_2/starting_ttl)  # proportion for murder

  if (progress == TRUE) {
    writeLines("\nBeginning sampling period.")
    sample_progress <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  }

  # Final Metropolis Hastings
  for (i in 1:iterations) {

    old_loglik <- fitMetrics(k_ends, full_data)  # calls fit matrix to have a function to start with

    k_and_q <- newEnds(k_ends, make_k, murder_k)
    k_ends_new <- k_and_q[[1]]
    q1 <- k_and_q[[2]]
    q2 <- k_and_q[[3]]
    type <- k_and_q[[4]]

    new_loglik <- fitMetrics(k_ends_new, full_data)

    # birth and death ratios

    product_runs <- c(0)

    old_new_product_birth_sampling <- function(product_runs) {
      for (j in k) {
        product_runs <- append(product_runs, (n - (3 * ar - (q1 + 1) * (2 *
          ar + j))))
      }
      return(prod(product_runs))
    }
    product_death <- function(k) {
      for (j in k) {
        product_runs <- append(product_runs, (n - 3 * ar - q1 * 2 * ar +
          j))
      }
      return(prod(product_runs))
    }
    old_new_product_death <- function(k) {
      for (j in k) {
        product_runs <- append(product_runs, (n - 3 * ar - (q1 - 1) * (2 *
          ar) + j))
      }
      return(prod(product_runs))
    }

    prior <- c()
    birth_old_to_new_ratio <- (factorial(starting_bkpts_2 + 1) * q2 * (dpois(length(k_ends_new) -
      2, lambda)))/(old_new_product_birth_sampling(starting_bkpts_2))
    birth_new_to_old_ratio <- (factorial(starting_bkpts_2) * q1 * (dpois(length(k_ends) -
      2, lambda)))/(product_death(starting_bkpts_2))
    death_new_to_old_ratio <- (factorial(starting_bkpts_2) * q2 * (dpois(length(k_ends) -
      2, lambda)))/(product_death(starting_bkpts_2))
    death_old_to_new_ratio <- (factorial(starting_bkpts_2 - 1) * make_k * (dpois(length(k_ends_new) -
      2, lambda)))/(old_new_product_death(starting_bkpts_2))
    ratios <- c(birth_old_to_new_ratio, birth_new_to_old_ratio, death_new_to_old_ratio,
      death_old_to_new_ratio)
    for (i in ratios) {
      # finds priors for birth and death ratios
      prior[i] = 1 - ratios[i]
    }
    # end of birth and death ratios

    delta_bic <- (-2 * new_loglik + log(n) * (length(k_ends_new) - 1) * (3 +
      ar)) - (-2 * old_loglik + log(n) * (length(k_ends) - 1) * (3 + ar))
    ratio <- (-1 * delta_bic/2) + (log(q1 * dpois(length(k_ends_new) - 2, lambda)) -
      log(q2 * dpois(length(k_ends) - 2, lambda)))
    u_ratio <- log(runif(1))  # random number from 0 to 1 taken from a uniform distribution and then log transformed

    if (abs(delta_bic) == Inf) {
      # safe guard against random models creating infinite ratios
      k_ends <- k_ends  # old
      bic <- (-2 * old_loglik + log(n) * (length(k_ends) - 1) * (3 + ar))
    } else if (ratio > u_ratio) {
      k_ends <- k_ends_new  # new
      bic <- (-2 * new_loglik + log(n) * (length(k_ends_new) - 1) * (3 + ar))
      accept_count <- accept_count + 1
      # looking at what type of step is done and accepted
      if (type == "add") {
        add.accept.count <<- add.accept.count + 1
      } else if (type == "sub") {
        sub.accept.count <<- sub.accept.count + 1
      } else if (type == "move") {
        move.accept.count <<- move.accept.count + 1
      } else if (type == "jiggle") {
        jiggle.accept.count <<- jiggle.accept.count + 1
      }
    } else {
      k_ends <- k_ends  # old
      bic <- (-2 * old_loglik + log(n) * (length(k_ends) - 1) * (3 + ar))
    }

    k <- k_ends[c(-1, -length(k_ends))]
    if (length(k) > ncol(all_k_best)) {
      all_k_best <- cbind(all_k_best, rep(NA, nrow(all_k_best)))
      all_k_best <- rbind(all_k_best, k)
    } else if (length(k) < ncol(all_k_best)) {
      k <- c(k, rep(NA, (ncol(all_k_best) - length(k))), recursive = T)
      all_k_best <- rbind(all_k_best, k)
    } else {
      all_k_best <- rbind(all_k_best, k)
    }
    all_BIC <- rbind(all_BIC, bic)

    if (fit_storage == TRUE) {
      fit_vec <- rep(NA_real_, n)
      squared_resids <- numeric(0)
      current_post_betas <- NULL
      current_post_sigmas <- NULL

      for (m in 2:length(k_ends)) {
        if (m > 2) {
          start_idx <- k_ends[m - 1] + 1
        } else {
          start_idx <- k_ends[m - 1]
        }
        end_idx <- k_ends[m]
        segment <- full_data[start_idx:end_idx, 2]
        if (length(segment) <= ar) {
          next
        }

        embed_mat <- stats::embed(segment, ar + 1)
        y_j <- embed_mat[, 1]
        lag_mat <- embed_mat[, -1, drop = FALSE]
        x_j <- cbind(1, lag_mat)

        ols_fit <- stats::lm.fit(x_j, y_j)
        sigma2_hat <- sum(ols_fit$residuals^2)/max(length(y_j) - ncol(x_j),
          1)
        if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
          sigma2_hat <- stats::var(y_j)
          if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
          sigma2_hat <- 1
          }
        }

        posterior_precision <- prior_precision + crossprod(x_j)/sigma2_hat
        v <- tryCatch(solve(posterior_precision), error = function(e) MASS::ginv(posterior_precision))
        beta_mean <- v %*% (prior_precision %*% b_0 + crossprod(x_j, y_j)/sigma2_hat)
        post_beta <- as.numeric(MASS::mvrnorm(1, as.numeric(beta_mean), v))

        fitted_vals <- as.numeric(x_j %*% post_beta)
        segment_fit <- c(rep(NA_real_, ar), fitted_vals)
        fit_vec[start_idx:end_idx] <- segment_fit
        squared_resids <- c(squared_resids, (y_j - fitted_vals)^2)

        resid_vec <- y_j - fitted_vals
        d0 <- as.numeric(0.5 * sum(resid_vec^2))
        rate <- ifelse(is.finite(d0) && d0 > 0, d0, .Machine$double.eps)
        shape <- (length(resid_vec)/2) + 2
        sigma_sample <- stats::rgamma(1, shape, rate = rate)
        post_sigma <- 1/sigma_sample

        post_beta_mat <- matrix(post_beta, ncol = 1)
        post_sigma_mat <- matrix(post_sigma, ncol = 1)
        current_post_betas <- if (is.null(current_post_betas))
          post_beta_mat else cbind(current_post_betas, post_beta_mat)
        current_post_sigmas <- if (is.null(current_post_sigmas))
          post_sigma_mat else cbind(current_post_sigmas, post_sigma_mat)
      }

      all_fits <- rbind(all_fits, fit_vec)
      all_MSE <- c(all_MSE, if (length(squared_resids)) mean(squared_resids) else NA_real_)

      beta_draws[[length(beta_draws) + 1]] <- if (is.null(current_post_betas)) {
        data.frame()
      } else {
        df <- as.data.frame(current_post_betas)
        rownames(df) <- paste0("B", 0:ar)
        colnames(df) <- paste0("Segment", seq_len(ncol(df)))
        df
      }

      sigma_draws[[length(sigma_draws) + 1]] <- if (is.null(current_post_sigmas)) {
        data.frame()
      } else {
        df <- as.data.frame(current_post_sigmas)
        rownames(df) <- "Sigma"
        colnames(df) <- paste0("Segment", seq_len(ncol(df)))
        df
      }
    }
  }

  if (fit_storage == TRUE) {
    if (length(all_MSE) > 0) {
      all_MSE <- data.frame(MSE = all_MSE)
    } else {
      all_MSE <- data.frame(MSE = numeric(0))
    }
    if (ncol(all_fits) > 0) {
      colnames(all_fits) <- c(1:ncol(all_fits))
      all_fits <- as.data.frame(all_fits)
    } else {
      all_fits <- data.frame()
    }
    post_beta_list <- beta_draws
    post_sigma_list <- sigma_draws
  }

  if (progress == TRUE) {
    writeLines("\n")
  }

  # cleaning up the matrices and counts
  if (length(all_k_best) != 0) {
    colnames(all_k_best) = c(1:ncol(all_k_best))
  }
  final.propose <- c(a.count, s.count, m.count, j.count)
  final.accept <- c(add.accept.count, sub.accept.count, move.accept.count, jiggle.accept.count)
  colnames(all_BIC) = "BIC"

  # getting distribution of k (number of breakpoints)
  if (nrow(all_k_best) > 0) {
    num_bkpts <- rowSums(!is.na(all_k_best))
  } else {
    num_bkpts <- numeric(0)
  }


  if (fit_storage == TRUE) {
    final_list <- list(accept_count/iterations, final.propose, final.accept,
      all_MSE, all_BIC, all_k_best, num_bkpts, post_beta_list, post_sigma_list,
      all_fits)
    names(final_list) = c("AcceptRate", "ProposedSteps", "AcceptedSteps", "MSE",
      "BIC", "Breakpoints", "NumBkpts", "Beta", "Sigma", "Fits")
  } else {
    final_list <- list(accept_count/iterations, final.propose, final.accept,
      all_BIC, all_k_best, num_bkpts)
    names(final_list) = c("AcceptRate", "ProposedSteps", "AcceptedSteps", "BIC",
      "Breakpoints", "NumBkpts")
  }

  return(final_list)
}


# calling the function test_data = test_data_44() current_result = baar(NA,
# test_data[,1], test_data[,2], 50, 50, jump=0.25, ar=1, progress=T,
# fit_storage=T)
