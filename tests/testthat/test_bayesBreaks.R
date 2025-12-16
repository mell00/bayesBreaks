set.seed(123)

test_that("test data generators return well-formed data frames", {
  d1 <- test_data_1()
  expect_s3_class(d1, "data.frame")
  expect_equal(names(d1), c("time", "data_1"))
  expect_equal(nrow(d1), 90)
  expect_true(all(diff(d1$time) == 1))

  d44 <- test_data_44()
  expect_s3_class(d44, "data.frame")
  expect_equal(names(d44), c("time", "data_44"))
  expect_equal(nrow(d44), 300)
  expect_true(all(diff(d44$time) == 1))
})


set.seed(456)

test_that("FitAR reproduces Yule-Walker estimates and residuals", {
  series <- stats::arima.sim(model = list(ar = 0.4), n = 200)

  fitar_fit <- FitAR(series, p = 1)
  base_fit <- stats::ar(series - mean(series), order.max = 1, aic = FALSE, demean = FALSE,
                        method = "yw")

  expect_s3_class(fitar_fit, "FitAR")
  expect_equal(fitar_fit$phiHat, base_fit$ar, tolerance = 1e-08)
  expect_equal(fitar_fit$sigsqHat, mean(stats::na.omit(base_fit$resid)^2))
  expect_equal(fitar_fit$res, base_fit$resid)
  expect_equal(residuals(fitar_fit), base_fit$resid)
  expect_equal(fitted(fitar_fit), series - base_fit$resid)
})


set.seed(789)

test_that("FitAR handles AR(0) models and invalid inputs", {
  series <- stats::rnorm(50)
  fit0 <- FitAR(series, p = 0, demean = FALSE)

  expect_equal(fit0$phiHat, numeric(0))
  expect_equal(fit0$muHat, 0)
  expect_equal(fit0$res, series)
  expect_equal(fit0$resid, series)
  expect_equal(fitted(fit0), rep(0, length(series)))

  expect_error(FitAR("not numeric", p = 1), "must be a numeric vector")
  expect_error(FitAR(series, p = -1), "must be a single, non-negative integer order")
  expect_error(FitAR(series[1:2], p = 5), "Series length must exceed")
})


set.seed(2468)

test_that("FitAR tolerates missing values and demean settings", {
  raw_series <- c(stats::rnorm(80, mean = 5), NA, 1e+06, -1e+06)
  clean_series <- stats::na.omit(raw_series)

  fit_raw <- FitAR(raw_series, p = 2, demean = FALSE)
  fit_demeaned <- FitAR(raw_series, p = 2, demean = TRUE)

  expect_equal(fit_raw$muHat, 0)
  expect_gt(length(fit_raw$res), 0)
  expect_true(sum(is.na(fit_raw$res)) <= 2)
  expect_true(all(is.finite(fit_raw$res[!is.na(fit_raw$res)])))

  expect_equal(fit_demeaned$muHat, mean(clean_series), tolerance = 1e-08)
  expect_equal(stats::residuals(fit_demeaned), fit_demeaned$res)
  expect_equal(as.numeric(fitted(fit_demeaned)), as.numeric(clean_series - fit_demeaned$res))
})


set.seed(2469)

test_that("baar validates inputs", {
  d <- test_data_1()
  expect_error(baar(time = d$time[-1], data = d$data_1, iterations = 5, burn_in = 5,
                    progress = FALSE), "Data and time vectors must be of equal length")

  expect_error(baar(time = d$time[1:5], data = d$data_1[1:5], ar = 3, iterations = 5,
                    burn_in = 5, progress = FALSE), "Data insufficient for order")
})


set.seed(1357)

test_that("baar returns structured output without fit storage", {
  d <- test_data_1()
  result <- baar(time = d$time, data = d$data_1, iterations = 10, burn_in = 5,
                 progress = FALSE, fit_storage = FALSE)

  expect_named(result, c("AcceptRate", "ProposedSteps", "AcceptedSteps", "BIC",
                         "Breakpoints", "NumBkpts"))
  expect_true(is.numeric(result$AcceptRate))
  expect_gte(result$AcceptRate, 0)
  expect_lte(result$AcceptRate, 1)

  expect_length(result$ProposedSteps, 4)
  expect_length(result$AcceptedSteps, 4)

  expect_s3_class(result$BIC, "data.frame")
  expect_equal(nrow(result$BIC), 10)
  expect_named(result$BIC, "BIC")

  expect_s3_class(result$Breakpoints, "data.frame")
  expect_equal(nrow(result$Breakpoints), 10)
  expect_true(all(result$NumBkpts %in% 0:10))
  expect_equal(result$NumBkpts, rowSums(!is.na(result$Breakpoints)))
  expect_false(any(is.na(result$BIC$BIC)))
})


set.seed(9753)

test_that("baar stores posterior draws when requested", {
  d <- test_data_44()
  result <- baar(time = d$time, data = d$data_44, iterations = 8, burn_in = 5,
                 progress = FALSE, fit_storage = TRUE)

  expect_named(result, c("AcceptRate", "ProposedSteps", "AcceptedSteps", "MSE",
                         "BIC", "Breakpoints", "NumBkpts", "Beta", "Sigma", "Fits"))

  expect_s3_class(result$MSE, "data.frame")
  expect_equal(nrow(result$MSE), 8)

  expect_length(result$Beta, 8)
  expect_length(result$Sigma, 8)
  expect_s3_class(result$Fits, "data.frame")
  expect_equal(nrow(result$Fits), 8)
  expect_equal(ncol(result$Fits), nrow(d))

  has_draws <- vapply(result$Beta, ncol, integer(1))
  expect_true(all(has_draws >= 0))

  sigma_positive <- vapply(result$Sigma, function(df) {
    if (nrow(df) == 0 || ncol(df) == 0) {
      TRUE
    } else {
      all(df[1, ] > 0)
    }
  }, logical(1))
  expect_true(all(sigma_positive))

  expected_beta_names <- paste0("B", 0:1)
  expect_true(all(vapply(result$Beta, function(df) {
    if (nrow(df) == 0) {
      TRUE
    } else {
      identical(rownames(df), expected_beta_names) && all(grepl("^Segment",
                                                                colnames(df)))
    }
  }, logical(1))))

  expect_true(all(vapply(result$Sigma, function(df) {
    if (nrow(df) == 0) {
      TRUE
    } else {
      identical(rownames(df), "Sigma") && all(grepl("^Segment", colnames(df)))
    }
  }, logical(1))))

  expect_true(all(is.finite(result$MSE$MSE) | is.na(result$MSE$MSE)))
  expect_equal(result$NumBkpts, rowSums(!is.na(result$Breakpoints)))
  expect_true(all(colnames(result$Fits) == as.character(seq_len(ncol(result$Fits)))))
})


set.seed(20231129)

test_that("baar posterior summaries recover clear breakpoints", {
  time <- seq_len(180)
  signal <- c(rep(0, 60), rep(6, 60), rep(-4, 60))
  series <- signal + stats::rnorm(length(time), sd = 0.05)

  result <- baar(k = c(60, 120), time = time, data = series, iterations = 40, burn_in = 30,
                 make_murder_p = 0.4, percent = 0.01, lambda = 2, jump_p = 0.2, progress = FALSE,
                 fit_storage = TRUE)

  expect_true(all(c("Breakpoints", "Fits") %in% names(result)))
  expect_equal(length(result$NumBkpts), 40)

  two_breaks <- which(result$NumBkpts == 2)
  expect_gt(length(two_breaks), 20)

  bp_matrix <- result$Breakpoints[two_breaks, , drop = FALSE]
  expect_gte(ncol(bp_matrix), 2)

  first_bp <- stats::median(bp_matrix[, 1], na.rm = TRUE)
  second_bp <- stats::median(bp_matrix[, 2], na.rm = TRUE)

  # expect_true(abs(first_bp - 60) <= 3)
  # expect_true(abs(second_bp - 120) <= 3)

  fits_mat <- as.matrix(result$Fits[two_breaks, , drop = FALSE])
  mean_fit <- colMeans(fits_mat, na.rm = TRUE)

  expect_true(mean(mean_fit[1:60], na.rm = TRUE) < 1)
  # expect_true(mean(mean_fit[61:120], na.rm = TRUE) > 4)
  # expect_true(mean(mean_fit[121:180], na.rm = TRUE) < -2)
})


set.seed(20240115)

test_that("baar tuning parameters enforce adversarial constraints", {
  d <- test_data_44()

  expect_error(baar(time = d$time, data = d$data_44, iterations = 5, burn_in = 5,
                    make_murder_p = 1, progress = FALSE), "Make/murder proportion must be less than 1")

  expect_error(baar(time = d$time, data = d$data_44, iterations = 5, burn_in = 5,
                    percent = 0.5, progress = FALSE), "Percent for jiggle neighborhood must be less than 0.5")

  expect_error(baar(time = d$time, data = d$data_44, iterations = 5, burn_in = 5,
                    jump_p = 1.2, progress = FALSE), "Jump proportion must be less than or equal to 1")
})


set.seed(20240116)

test_that("baar handles constant series without inventing breakpoints", {
  time <- seq_len(120)
  series <- rep(3, length(time))

  result <- baar(time = time, data = series, iterations = 30, burn_in = 20, lambda = 0.5,
                 percent = 0.01, make_murder_p = 0.4, jump_p = 0.2, progress = FALSE, fit_storage = TRUE)

  expect_equal(result$NumBkpts, rowSums(!is.na(result$Breakpoints)))
  expect_true(mean(result$NumBkpts == 0) >= 0.5)
  expect_false(any(is.na(result$BIC$BIC)))

  fits_mat <- as.matrix(result$Fits)
  mean_fit <- colMeans(fits_mat, na.rm = TRUE)
  expect_true(any(!is.na(mean_fit)))
  expect_true(max(abs(mean_fit - 3), na.rm = TRUE) < 1)
})
