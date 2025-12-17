#' Rcpp-accelerated BAAR sampler
#'
#' C++ implementation of the Bayesian Adaptive Auto-Regression sampler,
#' mirroring the interface of [baar()] while relying on Armadillo routines for
#' likelihood evaluation and proposals.
#'
#' @inheritParams baar
#' @export
baar_rcpp <- function(k = NULL, time, data, iterations, burn_in = 50,
                      make_murder_p = 0.5, percent = 0.02, lambda = 1, jump_p = 0.25,
                      ar = 1, progress = TRUE, fit_storage = FALSE) {
  out <- .Call("rcpp_baar", PACKAGE = "bayesBreaks", k, time, data, as.integer(iterations),
               as.integer(burn_in), as.numeric(make_murder_p), as.numeric(percent),
               as.numeric(lambda), as.numeric(jump_p), as.integer(ar), as.logical(progress),
               as.logical(fit_storage))

  # similar to R implementation w/ data.frame outputs and numeric BIC column
  out$BIC <- data.frame(BIC = as.numeric(out$BIC))

  bp_mat <- as.matrix(out$Breakpoints)
  if (ncol(bp_mat) == 0) {
    out$Breakpoints <- data.frame(matrix(nrow = iterations, ncol = 0))
  } else {
    colnames(bp_mat) <- as.character(seq_len(ncol(bp_mat)))
    out$Breakpoints <- as.data.frame(bp_mat)
  }

  out
}

