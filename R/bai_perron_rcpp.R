#' Rcpp Bai-Perron breakpoint search
#'
#' C++ implementation of the Bai-Perron dynamic-programming routine
#' used by [bai_perron_ar()], returning the same breakpoint summary but using
#' compiled Armadillo linear algebra for performance.
#'
#' @inheritParams bai_perron_ar
#' @return List containing breakpoints, all_breakpoints, SSR, and BIC vectors
#'   matching the native R implementation.
#' @export
bai_perron_ar_rcpp <- function(data, order = 1L, interval = 0.15, max_breaks = 3L) {
  .Call("rcpp_bai_perron", PACKAGE = "bayesBreaks", data, as.integer(order),
        as.numeric(interval), as.integer(max_breaks))
}
