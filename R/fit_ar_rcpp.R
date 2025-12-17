#' Autoregressive fitting via RcppArmadillo
#'
#' Alternative to [FitAR()] using Yule-Walker equations and
#' Armadillo linear algebra.
#'
#' @inheritParams FitAR
#' @export
FitAR_rcpp <- function(z, p, demean = TRUE) {
  res <- .Call("rcpp_fit_ar", PACKAGE = "bayesBreaks", z, as.integer(p), as.logical(demean))
  tsp_in <- stats::tsp(z)
  start_val <- if (!is.null(tsp_in)) tsp_in[1] else 1
  freq_val <- if (!is.null(tsp_in)) tsp_in[3] else 1

  res$res <- stats::ts(as.numeric(res$res), start = start_val, frequency = freq_val)
  res$fits <- stats::ts(as.numeric(res$fits), start = start_val, frequency = freq_val)

  # drop matrix attributes for parity with FitAR()
  res$phiHat <- as.numeric(res$phiHat)


  res
}

