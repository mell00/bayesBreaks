#' Autoregressive fitting via RcppArmadillo
#'
#' Alternative to [FitAR()] using Yule-Walker equations and
#' Armadillo linear algebra.
#'
#' @inheritParams FitAR
#' @export
FitAR_rcpp <- function(z, p, demean = TRUE) {
  .Call("rcpp_fit_ar", PACKAGE = "bayesBreaks", z, as.integer(p), as.logical(demean))
}
