

simulate_segments <- function(means, sds, lengths, column_name) {
  values <- unlist(Map(function(mean, sd, len) stats::rnorm(len, mean = mean, sd = sd),
                       means, sds, lengths))
  time <- seq_along(values)
  data.frame(time = time, stats::setNames(list(values), column_name), check.names = FALSE)
}

test_data_0_a <- function() {
  simulate_segments(10, 1, 90, "data_0_a")
}

test_data_0_b <- function() {
  simulate_segments(10, 5, 90, "data_0_b")
}

test_data_1 <- function() {
  simulate_segments(c(10, 20), c(1, 1), c(45, 45), "data_1")
}

test_data_2 <- function() {
  simulate_segments(c(5, 15, 30), c(1, 1, 1), c(30, 30, 30), "data_2")
}

test_data_3 <- function() {
  simulate_segments(c(10, 20), c(5, 5), c(45, 45), "data_3")
}

test_data_4 <- function() {
  simulate_segments(c(5, 20, 35), c(6, 6, 6), c(30, 30, 30), "data_4")
}

test_data_5 <- function() {
  first <- stats::rnorm(45, mean = seq_len(45), sd = 1)
  second <- stats::rnorm(45, mean = seq(46, 224, by = 4), sd = 1)
  data <- c(first, second)
  time <- seq_along(data)
  data.frame(time = time, data_5 = data)
}

test_data_6 <- function() {
  first <- stats::rnorm(30, mean = seq_len(30), sd = 1)
  second <- stats::rnorm(30, mean = seq(31, 148, by = 4), sd = 1)
  third <- stats::rnorm(30, mean = seq(149, 382, by = 8), sd = 1)
  data <- c(first, second, third)
  time <- seq_along(data)
  data.frame(time = time, data_6 = data)
}

test_data_7 <- function() {
  first <- stats::rnorm(length(seq(1, 23, by = 0.5)), mean = seq(1, 23, by = 0.5),
                        sd = 1)
  second <- stats::rnorm(length(24:68), mean = 24:68, sd = 1)
  data <- c(first, second)
  time <- seq_along(data)
  data.frame(time = time, data_7 = data)
}

test_data_8 <- function() {
  first_seq <- seq(1, 8.25, by = 0.25)
  second_seq <- seq(8.5, 23, by = 0.5)
  third_seq <- seq(24, 46, by = 0.75)
  first <- stats::rnorm(length(first_seq), mean = first_seq, sd = 1)
  second <- stats::rnorm(length(second_seq), mean = second_seq, sd = 1)
  third <- stats::rnorm(length(third_seq), mean = third_seq, sd = 1)
  data <- c(first, second, third)
  time <- seq_along(data)
  data.frame(time = time, data_8 = data)
}

test_data_9 <- function() {
  simulate_segments(c(10, 10), c(1, 5), c(45, 45), "data_9")
}

test_data_10 <- function() {
  simulate_segments(c(10, 10, 10), c(1, 5, 1), c(30, 30, 30), "data_10")
}

test_data_11 <- function() {
  first <- stats::arima.sim(model = list(ar = 0.1), n = 45)
  second <- stats::arima.sim(model = list(ar = 0.9999), n = 45)
  data <- c(second, first)
  time <- seq_along(data)
  data.frame(time = time, data_11 = data)
}

test_data_44 <- function() {
  simulate_segments(c(10, 20, 30), c(5, 5, 5), c(100, 100, 100), "data_44")
}

test_data_100 <- function() {
  segments <- rep(c(5, 10, 15, 10, 5, 10, 15, 10, 5), each = 1)
  means <- segments
  sds <- rep(1, length(means))
  lengths <- rep(100, length(means))
  simulate_segments(means, sds, lengths, "data_100")
}

test_data_200 <- function() {
  means <- rep(10, 9)
  sds <- c(10, 5, 1, 10, 1, 10, 10, 5, 1)
  lengths <- rep(100, length(means))
  simulate_segments(means, sds, lengths, "data_200")
}

test_data_300 <- function() {
  simulate_segments(c(10, 20), c(5, 5), c(100, 100), "data_300")
}
