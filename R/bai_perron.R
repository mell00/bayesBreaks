
#' Bai-Perron breakpoint search for autoregressive segments
#'
#' This helper mirrors the classical Bai-Perron dynamic-programming routine
#' that shipped with the REU BAAR scripts.  It evaluates all admissible
#' segmentations of an autoregressive time series up to a supplied maximum
#' number of breakpoints, fits an AR model within each segment using the
#' lightweight [`FitAR()`][FitAR], and selects the configuration that minimises
#' the Bayesian information criterion (BIC).
#'
#' @param data Numeric vector containing the observed time-series values.
#' @param order Non-negative integer giving the autoregressive order to fit in
#'   each segment.
#' @param interval Minimum proportion of observations that must fall inside any
#'   segment.  This mirrors the original Bai-Perron `epsilon` constraint and
#'   prevents degenerate breakpoint placement.
#' @param max_breaks Maximum number of breakpoints to consider.
#' @param progress Logical; if `TRUE` a text progress bar summarises the
#'   pre-computation of segment residual sums of squares.
#'
#' @return A list containing
#'   \describe{
#'     \item{breakpoints}{The breakpoint locations (indices) for the model with
#'       the smallest BIC.  An empty integer vector indicates that the null model
#'       with zero breaks was preferred.}
#'     \item{all_breakpoints}{A list where the `k`-th element stores the
#'       breakpoint locations associated with the optimal `k`-break model.}
#'     \item{SSR}{A named numeric vector giving the residual sum of squares for
#'       the optimal configuration with `k` breakpoints, including the `0`
#'       breakpoint null model.}
#'     \item{BIC}{A named numeric vector giving the Bayesian information
#'       criterion for each breakpoint count.}
#'   }
#'
#' @examples
#' set.seed(1)
#' series <- test_data_2()
#' bai_perron_ar(series$data_2, order = 0, interval = 0.1, max_breaks = 2,
#'   progress = FALSE)
#'
#' @export

bai_perron_ar <- function(data, order = 1L, interval = 0.15, max_breaks = 3L,
                          progress = interactive()) {
  if (!is.numeric(data)) {
    stop("`data` must be a numeric vector.")
  }
  if (length(order) != 1L || order < 0 || floor(order) != ceiling(order)) {
    stop("`order` must be a single non-negative integer.")
  }
  if (!is.numeric(interval) || length(interval) != 1L || interval <= 0 ||
      interval >= 1) {
    stop("`interval` must be a single number between 0 and 1.")
  }
  if (length(max_breaks) != 1L || max_breaks < 0 ||
      floor(max_breaks) != ceiling(max_breaks)) {
    stop("`max_breaks` must be a single non-negative integer.")
  }
  if (!is.logical(progress) || length(progress) != 1L) {
    stop("`progress` must be a single logical value.")
  }

  y <- stats::na.omit(as.numeric(data))
  n <- length(y)
  if (n <= order) {
    stop("Series length must exceed the autoregressive order.")
  }

  min_segment <- max(order + 1L, floor(n * interval))
  min_segment <- max(min_segment, 3L)
  if (min_segment > n) {
    stop("`interval` is too restrictive for the provided series length.")
  }

  if (min_segment * (max_breaks + 1L) > n) {
    max_allowed <- max(floor(n / min_segment) - 1L, 0L)
    stop(sprintf("`max_breaks` is too large for the requested interval.  Try %d or fewer.",
      max_allowed))
  }

  compute_ssr <- function(start, end) {
    if ((end - start + 1L) <= order) {
      return(Inf)
    }
    fit <- FitAR(y[start:end], p = order)
    res <- stats::na.omit(residuals(fit))
    sum(res^2)
  }

  ssr_matrix <- matrix(Inf, nrow = n, ncol = n)
  idx <- seq_len(n - min_segment + 1L)
  pb <- NULL
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = length(idx), style = 3)
  }
  for (i in idx) {
    start <- i
    end_seq <- seq.int(start + min_segment - 1L, n)
    for (end in end_seq) {
      ssr_matrix[start, end] <- compute_ssr(start, end)
    }
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (!is.null(pb)) {
    close(pb)
  }

  max_segments <- max_breaks + 1L
  dp <- matrix(Inf, nrow = max_segments, ncol = n)
  last_break <- matrix(NA_integer_, nrow = max_segments, ncol = n)

  for (end in seq.int(min_segment, n)) {
    ssr <- ssr_matrix[1, end]
    if (is.finite(ssr)) {
      dp[1, end] <- ssr
      last_break[1, end] <- 0L
    }
  }

  if (max_segments > 1L) {
    for (segments in 2:max_segments) {
      min_end <- segments * min_segment
      if (min_end > n) {
        break
      }
      for (end in seq.int(min_end, n)) {
        best_total <- Inf
        best_break <- NA_integer_
        candidate_breaks <- seq.int((segments - 1L) * min_segment,
          end - min_segment)
        for (break_point in candidate_breaks) {
          prev_total <- dp[segments - 1L, break_point]
          if (!is.finite(prev_total)) {
            next
          }
          seg_ssr <- ssr_matrix[break_point + 1L, end]
          if (!is.finite(seg_ssr)) {
            next
          }
          total <- prev_total + seg_ssr
          if (total < best_total) {
            best_total <- total
            best_break <- break_point
          }
        }
        if (is.finite(best_total)) {
          dp[segments, end] <- best_total
          last_break[segments, end] <- best_break
        }
      }
    }
  }

  reconstruct_breaks <- function(segments) {
    if (!is.finite(dp[segments, n])) {
      return(integer(0))
    }
    breaks <- integer(segments - 1L)
    current_end <- n
    current_segments <- segments
    while (current_segments > 1L) {
      bp <- last_break[current_segments, current_end]
      if (is.na(bp) || bp <= 0L) {
        return(integer(0))
      }
      breaks[current_segments - 1L] <- bp
      current_end <- bp
      current_segments <- current_segments - 1L
    }
    breaks
  }

  ssr_values <- numeric(max_breaks + 1L)
  bic_values <- numeric(max_breaks + 1L)
  breakpoint_sets <- vector("list", length = max_breaks)

  null_fit <- FitAR(y, p = order)
  null_res <- stats::na.omit(residuals(null_fit))
  null_ssr <- sum(null_res^2)
  base_constant <- n * (log(2 * pi) + 1)
  ssr_values[1] <- null_ssr
  bic_values[1] <- n * log(null_ssr / n) + base_constant +
    log(n) * (order + 2) * 1L

  for (breaks in seq_len(max_breaks)) {
    segments <- breaks + 1L
    ssr_total <- dp[segments, n]
    if (!is.finite(ssr_total)) {
      ssr_values[breaks + 1L] <- NA_real_
      bic_values[breaks + 1L] <- NA_real_
      breakpoint_sets[[breaks]] <- integer(0)
      next
    }
    ssr_values[breaks + 1L] <- ssr_total
    bic_values[breaks + 1L] <- n * log(ssr_total / n) + base_constant +
      log(n) * (order + 2) * segments
    breakpoint_sets[[breaks]] <- reconstruct_breaks(segments)
  }

  names(ssr_values) <- as.character(0:max_breaks)
  names(bic_values) <- as.character(0:max_breaks)
  names(breakpoint_sets) <- as.character(seq_len(max_breaks))

  best_index <- which.min(bic_values)
  best_breakpoints <- if (best_index == 1L) {
    integer(0)
  } else {
    breakpoint_sets[[best_index - 1L]]
  }

  list(
    breakpoints = best_breakpoints,
    all_breakpoints = breakpoint_sets,
    SSR = ssr_values,
    BIC = bic_values
  )
}
