#' Package initialization
#'
#' Registers package's compiled code. imports helpers needed for routine usage.
#' @useDynLib bayesBreaks, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
