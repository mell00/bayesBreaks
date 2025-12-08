
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
  # Validate basic input types and shapes
  if (!is.numeric(z)) {
    stop("`z` must be a numeric vector.")
  }
  if (length(p) != 1L || p < 0 || floor(p) != ceiling(p)) {
    stop("`p` must be a single, non-negative integer order.")
  }

  tsp_z <- stats::tsp(z)
  with_tsp <- function(x) {
    if (is.null(tsp_z) || length(x) == 0L) {
      return(x)
    }
    stats::ts(as.numeric(x), start = tsp_z[1], frequency = tsp_z[3])
  }

  # Drop missing values and coerce to a numeric vector
  z <- as.numeric(z)
  z <- stats::na.omit(z)
  n <- length(z)
  if (n <= p) {
    stop("Series length must exceed the autoregressive order.")
  }


  # (optional) remove the sample mean before fitting AR coefficients
  mu_hat <- if (demean) {
    mean(z)
  } else {
    0
  }
  y <- z - mu_hat

  # handle AR(0) case directly to avoid unnecessary model fitting
  if (p == 0L) {
    res <- with_tsp(y)
    sigma2 <- if (n > 0) {
      mean(res^2)
    } else {
      NA_real_
    }
    if (is.na(sigma2) || sigma2 <= 0) {
      sigma2 <- .Machine$double.eps
    }
    loglik <- -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
    fits <- with_tsp(z - res)
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

  # short-circuit constant series to avoid Yule-Walker failures
  if (stats::var(y) == 0) {
    phi <- rep(0, p)
    res <- with_tsp(rep(0, n))
    sigma2 <- .Machine$double.eps
    loglik <- -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
    fits <- with_tsp(z)
  } else {
    # estimate AR(p) coefficients using Yule-Walker equations
    fit <- stats::ar(y, order.max = p, aic = FALSE, demean = FALSE, method = "yw")
    phi <- fit$ar
    if (length(phi) < p) {
      phi <- c(phi, rep(0, p - length(phi)))
    }

    # summarize residual scale and implied Gaussian log-likelihood
    res <- with_tsp(fit$resid)
    sigma2 <- mean(stats::na.omit(as.numeric(res))^2)
    if (is.na(sigma2) || sigma2 <= 0) {
      sigma2 <- .Machine$double.eps
    }
    loglik <- -0.5 * n * (log(2 * pi) + log(sigma2) + 1)
    fits <- with_tsp(z - res)
  }

  ans <- list(loglikelihood = loglik, phiHat = phi, sigsqHat = sigma2, muHat = mu_hat,
              covHat = matrix(numeric(0), nrow = 0, ncol = 0), zetaHat = phi, RacfMatrix = matrix(numeric(0),
                                                                                                  nrow = 0, ncol = 0), LjungBoxQ
              = NULL, res = res, resid = res, fits = fits,
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
