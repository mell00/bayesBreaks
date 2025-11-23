
#' FitAR package
#'
#' The original BAAR scripts relied on the `FitAR::FitAR()` function to fit
#' autoregressive models within each breakpoint segment. The CRAN package was
#' archived, so this version re-implements the subset
#' of functionality required by those scripts. It fits a standard
#' autoregressive model of the requested order using Yule-Walker equations and
#' returns an object that mimics the structure of the historical FitAR output.
#'
#' @param z Numeric vector containing the observed time-series values.
#' @param p Positive integer giving the autoregressive order. Only scalar orders
#'   are supported by this compatibility shim.
#' @param demean Logical; if `TRUE` the series mean is removed prior to
#'   estimating the AR coefficients.
#' @param ... Ignored, included for API compatibility with the retired
#'   `FitAR` package.
#'
#' @return A list with class `'FitAR'` containing components compatible with
#'   the original package, including the fitted coefficients, residuals, and the
#'   Gaussian log-likelihood under the estimated innovation variance.
#'
#' @examples
#' FitAR(stats::arima.sim(model = list(ar = 0.5), n = 100), 1)
#'
#' @export
#'
FitAR <- function(z, p, demean = TRUE, ...) {
  if (!is.numeric(z)) {
    stop("`z` must be a numeric vector.")
  }
  if (length(p) != 1L || p < 0 || floor(p) != ceiling(p)) {
    stop("`p` must be a single, non-negative integer order.")
  }

  z <- as.numeric(z)
  z <- stats::na.omit(z)
  n <- length(z)
  if (n <= p) {
    stop("Series length must exceed the autoregressive order.")
  }

  mu_hat <- if (demean) {
    mean(z)
  } else {
    0
  }
  y <- z - mu_hat

  if (p == 0L) {
    res <- y
    sigma2 <- if (n > 0) {
      mean(res^2)
    } else {
      NA_real_
    }
    loglik <- if (!is.na(sigma2) && sigma2 > 0) {
      -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
    } else {
      NA_real_
    }
    fits <- z - res
    ans <- list(loglikelihood = loglik, phiHat = numeric(0), sigsqHat = sigma2,
                muHat = mu_hat, covHat = matrix(numeric(0), nrow = 0, ncol = 0), zetaHat = numeric(0),
                RacfMatrix = matrix(numeric(0), nrow = 0, ncol = 0), LjungBoxQ = NULL,
                res = res, resid = res, fits = fits, SubsetQ = FALSE, pvec = integer(0),
                demean = demean, FitMethod = "YW", iterationCount = 0L, convergence = 0L,
                MeanMLE = FALSE, tsp = stats::tsp(z), call = match.call(), ARModel = "ARz",
                DataTitle = attr(z, "title"), ModelTitle = "AR(0)", z = z)
    class(ans) <- "FitAR"
    return(ans)
  }

  fit <- stats::ar(y, order.max = p, aic = FALSE, demean = FALSE, method = "yw")
  phi <- fit$ar
  if (length(phi) < p) {
    phi <- c(phi, rep(0, p - length(phi)))
  }

  res <- fit$resid
  sigma2 <- mean(stats::na.omit(res)^2)
  loglik <- -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
  fits <- z - res

  ans <- list(loglikelihood = loglik, phiHat = phi, sigsqHat = sigma2, muHat = mu_hat,
              covHat = matrix(numeric(0), nrow = 0, ncol = 0), zetaHat = phi, RacfMatrix = matrix(numeric(0),
                                                                                                  nrow = 0, ncol = 0), LjungBoxQ = NULL, res = res, resid = res, fits = fits,
              SubsetQ = FALSE, pvec = seq_len(p), demean = demean, FitMethod = "YW", iterationCount = 0L,
              convergence = 0L, MeanMLE = FALSE, tsp = stats::tsp(z), call = match.call(),
              ARModel = "ARz", DataTitle = attr(z, "title"), ModelTitle = sprintf("AR(%d)",
                                                                                  p), z = z)
  class(ans) <- "FitAR"
  ans
}

#' @export
residuals.FitAR <- function(object, ...) {
  object$res
}

#' @export
fitted.FitAR <- function(object, ...) {
  object$fits
}
