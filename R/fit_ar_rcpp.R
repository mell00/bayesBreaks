#' Autoregressive fitting via RcppArmadillo
#'
#' Alternative to [FitAR()] using Yule-Walker equations and
#' Armadillo linear algebra.
#'
#' @inheritParams FitAR
#' @export
FitAR_rcpp <- function(z, p, demean = TRUE) {
  res <- .Call("rcpp_fit_ar", PACKAGE = "bayesBreaks", z, as.integer(p), as.logical(demean))

  # drop any matrix dimensions introduced by Armadillo and restore ts attributes
  res$phiHat <- as.numeric(res$phiHat)
  res$res <- as.numeric(res$res)
  res$fits <- as.numeric(res$fits)

  tsp_in <- stats::tsp(z)
  if (!is.null(tsp_in)) {
    res$res <- stats::ts(res$res, start = tsp_in[1], frequency = tsp_in[3])
    res$fits <- stats::ts(res$fits, start = tsp_in[1], frequency = tsp_in[3])
  }

  res
}

