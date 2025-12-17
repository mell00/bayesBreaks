set.seed(20231201)

test_that("Bai-Perron helper mirrors classical search", {
  data <- test_data_2()
  result <- bai_perron_ar(data$data_2, order = 0, interval = 0.1, max_breaks = 2,
                          progress = FALSE)

  expect_named(result, c("breakpoints", "all_breakpoints", "SSR", "BIC"))

  expect_length(result$breakpoints, 2)
  expect_true(all(abs(result$breakpoints - c(30, 60)) <= 5))

  expect_equal(length(result$all_breakpoints), 2)
  expect_true(all(vapply(result$all_breakpoints, is.integer, logical(1))))

  expect_length(result$SSR, 3)
  expect_named(result$SSR, c("0", "1", "2"))

  expect_length(result$BIC, 3)
  expect_named(result$BIC, c("0", "1", "2"))

  expect_true(which.min(result$BIC) >= 2)
})


test_that("Bai-Perron helper validates adversarial inputs", {
  expect_error(bai_perron_ar("not numeric"), "must be a numeric vector")
  expect_error(bai_perron_ar(1:10, order = -1), "non-negative integer")
  expect_error(bai_perron_ar(1:10, interval = 0), "between 0 and 1")
  expect_error(bai_perron_ar(1:10, max_breaks = -2), "non-negative integer")
  expect_error(bai_perron_ar(1:5, order = 3), "Series length must exceed")
  expect_error(bai_perron_ar(rep(0, 50), interval = 0.9, max_breaks = 3),
               "too restrictive")
})

test_that("Bai-Perron supports alternative grids and AR orders", {
  data <- test_data_2()

  fine_grid <- bai_perron_ar(data$data_2, order = 0, interval = 0.05, max_breaks = 2,
                             progress = FALSE)
  expect_length(fine_grid$breakpoints, 2)
  expect_true(all(abs(fine_grid$breakpoints - c(30, 60)) <= 6))

  ar_signal <- arima.sim(model = list(ar = 0.5), n = 160)
  ar_series <- ar_signal + rep(c(0, 3), each = 80)
  ar_result <- bai_perron_ar(ar_series, order = 1, interval = 0.1, max_breaks = 1,
                             progress = FALSE)

  expect_length(ar_result$breakpoints, 1)
  expect_true(abs(ar_result$breakpoints - 80) <= 6)
  expect_true(which.min(ar_result$BIC) >= 1)
})



test_that("Rcpp Bai-Perron mirrors R results and runtime", {
  skip_on_cran()

  data <- test_data_2()

  r_out <- bai_perron_ar(data$data_2, order = 0, interval = 0.12, max_breaks = 2,
                         progress = FALSE)
  c_out <- bai_perron_ar_rcpp(data$data_2, order = 0, interval = 0.12, max_breaks = 2)

  expect_equal(c_out$breakpoints, r_out$breakpoints)
  expect_equal(c_out$all_breakpoints, r_out$all_breakpoints)
  expect_equal(c_out$SSR, r_out$SSR, tolerance = 1e-08)
  expect_equal(c_out$BIC, r_out$BIC, tolerance = 1e-08)

  r_time <- system.time(bai_perron_ar(data$data_2, order = 0, interval = 0.12, max_breaks = 2,
                                      progress = FALSE))["elapsed"]
  c_time <- system.time(bai_perron_ar_rcpp(data$data_2, order = 0, interval = 0.12, max_breaks = 2))["elapsed"]

  expect_lte(c_time, r_time * 2)
})
