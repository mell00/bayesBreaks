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
  res <- .Call("rcpp_bai_perron", PACKAGE = "bayesBreaks", data, as.integer(order),
               as.numeric(interval), as.integer(max_breaks))

  # normalise names and fall back to the R helper if the compiled result is degenerate
  names(res$SSR) <- names(res$BIC) <- as.character(0:max_breaks)

  bp <- res$breakpoints
  all_bp <- res$all_breakpoints
  if ((max_breaks > 0 && length(bp) > 0 && all(bp == 0)) ||
      length(bp) == 0 ||
      (max_breaks > 0 && length(all_bp) > 0 && all(vapply(all_bp, length, integer(1)) == 0))) {
    return(bai_perron_ar(data, order = order, interval = interval, max_breaks = max_breaks, progress = FALSE))
  }

  res
}
