#' Dummy compiled helper
#'
#' Provides a minimal compiled algorithm backed by Rcpp so the package keeps
#' working native dependency. Returns a constant value from
#' compiled code path.
#'
#' @return Integer scalar (always 42).
#' @export
rcpp_dummy_value <- function() {
  .Call("rcpp_dummy_value", PACKAGE = "bayesBreaks")
}
