#' Simulated BAAR benchmarking datasets
#'
#' Collection of synthetic time-series data frames used to validate and
#' benchmark the Bayesian Adaptive Auto-Regression (BAAR) sampler. Each dataset
#' generated from the helper functions documented in
#' \code{\link{test_data_generators}} and stored for reuse in
#' examples/vignettes. Every object contains a `time` column of sequential
#' indices and a second column named after the dataset (for example
#' `data_0_a`).
#'
#' @format A data frame with two columns:
#' \describe{
#'   \item{time}{Integer sequence denoting the observation index.}
#'   \item{<dataset>}{Numeric observations for the simulated series. The column
#'   name matches the dataset object (e.g., `data_10`).}
#' }
#' @source Generated via functions in \code{\link{test_data_generators}}.
#' @aliases data_0_a data_0_b data_1 data_2 data_3 data_4 data_5 data_6 data_7
#' data_8 data_9 data_10 data_11 data_44 data_100 data_200 data_300
#' @docType data
#' @keywords datasets
#' @examples
#' data(data_0_a)
#' head(data_0_a)
"data_0_a"

#' @rdname data_0_a
"data_0_b"

#' @rdname data_0_a
"data_1"

#' @rdname data_0_a
"data_2"

#' @rdname data_0_a
"data_3"

#' @rdname data_0_a
"data_4"

#' @rdname data_0_a
"data_5"

#' @rdname data_0_a
"data_6"

#' @rdname data_0_a
"data_7"

#' @rdname data_0_a
"data_8"

#' @rdname data_0_a
"data_9"

#' @rdname data_0_a
"data_10"

#' @rdname data_0_a
"data_11"

#' @rdname data_0_a
"data_44"

#' @rdname data_0_a
"data_100"

#' @rdname data_0_a
"data_200"

#' @rdname data_0_a
"data_300"
