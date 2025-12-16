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
