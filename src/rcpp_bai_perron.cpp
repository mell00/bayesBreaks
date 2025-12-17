#include <RcppArmadillo.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::vec;
using arma::uword;

// simple AR fit returning SSR
static double segment_ssr(const vec &y, int order) {
  if (order == 0) {
    return arma::accu(arma::square(y - arma::mean(y)));
  }
  uword n = y.n_elem;
  if (n <= static_cast<uword>(order)) return std::numeric_limits<double>::infinity();

  double mean_y = arma::mean(y);
  vec centered = y - mean_y;
  // autocovariances
  vec gamma(order + 1, arma::fill::zeros);
  for (int lag = 0; lag <= order; ++lag) {
    double acc = 0.0;
    for (uword i = 0; i + lag < n; ++i) acc += centered[i] * centered[i + lag];
    gamma[lag] = acc / static_cast<double>(n);
  }
  arma::mat toeplitz(order, order);
  for (int i = 0; i < order; ++i) {
    for (int j = 0; j < order; ++j) toeplitz(i, j) = gamma[std::abs(i - j)];
  }
  vec rhs = gamma.tail(order);
  vec phi = arma::solve(toeplitz, rhs);
  vec resid(n, arma::fill::zeros);
  for (uword t = 0; t < n; ++t) {
    double pred = 0.0;
    for (int j = 0; j < order; ++j) {
      if (t <= static_cast<uword>(j)) break;
      pred += phi[j] * centered[t - j - 1];
    }
    resid[t] = centered[t] - pred;
  }
  return arma::accu(arma::square(resid));
}

extern "C" SEXP rcpp_bai_perron(SEXP data, SEXP order, SEXP interval, SEXP max_breaks) {
  vec y = Rcpp::as<vec>(data);
  y = y.elem(arma::find_finite(y));
  uword n = y.n_elem;
  if (n == 0) Rcpp::stop("`data` must contain observations.");
  int order_val = Rcpp::as<int>(order);
  double interval_val = Rcpp::as<double>(interval);
  int max_breaks_val = Rcpp::as<int>(max_breaks);
  if (order_val < 0) Rcpp::stop("`order` must be non-negative.");
  if (interval_val <= 0 || interval_val >= 1) Rcpp::stop("`interval` must be in (0, 1).");
  if (max_breaks_val < 0) Rcpp::stop("`max_breaks` must be non-negative.");

  uword min_segment = std::max<std::size_t>(order_val + 1, std::floor(n * interval_val));
  min_segment = std::max<uword>(min_segment, 3);
  if (min_segment > n || min_segment * static_cast<uword>(max_breaks_val + 1) > n) {
    Rcpp::stop("`interval` is too restrictive for the provided series length.");
  }

  arma::mat ssr_matrix(n, n, arma::fill::value(std::numeric_limits<double>::infinity()));
  for (uword start = 0; start < n; ++start) {
    uword min_end = start + min_segment - 1;
    if (min_end >= n) break;
    for (uword end = min_end; end < n; ++end) {
      ssr_matrix(start, end) = segment_ssr(y.subvec(start, end), order_val);
    }
  }

  int max_segments = max_breaks_val + 1;
  arma::mat dp(max_segments, n, arma::fill::value(std::numeric_limits<double>::infinity()));
  arma::imat last_break(max_segments, n, arma::fill::value(-1));

  for (uword end = min_segment - 1; end < n; ++end) {
    double ssr = ssr_matrix(0, end);
    if (R_finite(ssr)) {
      dp(0, end) = ssr;
      last_break(0, end) = 0;
    }
  }

  for (int seg = 1; seg < max_segments; ++seg) {
    uword min_end = (seg + 1) * min_segment - 1;
    if (min_end >= n) break;
    for (uword end = min_end; end < n; ++end) {
      double best_total = std::numeric_limits<double>::infinity();
      int best_break = -1;
      uword start_min = seg * min_segment - 1;
      uword start_max = end - min_segment;
      for (uword br = start_min; br <= start_max; ++br) {
        double prev = dp(seg - 1, br);
        if (!R_finite(prev)) continue;
        double seg_ssr = ssr_matrix(br + 1, end);
        if (!R_finite(seg_ssr)) continue;
        double total = prev + seg_ssr;
        if (total < best_total) {
          best_total = total;
          best_break = br + 1; // store start of segment as breakpoint index
        }
      }
      if (R_finite(best_total)) {
        dp(seg, end) = best_total;
        last_break(seg, end) = best_break;
      }
    }
  }

  Rcpp::NumericVector ssr_values(max_breaks_val + 1, NA_REAL);
  Rcpp::NumericVector bic_values(max_breaks_val + 1, NA_REAL);
  std::vector< Rcpp::IntegerVector > breakpoint_sets(max_breaks_val);
  double base_constant = static_cast<double>(n) * (std::log(2 * M_PI) + 1.0);

  // null model
  double null_ssr = segment_ssr(y, order_val);
  ssr_values[0] = null_ssr;
  bic_values[0] = static_cast<double>(n) * std::log(null_ssr / static_cast<double>(n)) +
    base_constant + std::log(static_cast<double>(n)) * (order_val + 2);

  for (int breaks = 1; breaks <= max_breaks_val; ++breaks) {
    int segments = breaks + 1;
    double ssr_total = dp(segments - 1, n - 1);
    if (!R_finite(ssr_total)) {
      breakpoint_sets[breaks - 1] = Rcpp::IntegerVector();
      continue;
    }
    Rcpp::IntegerVector bps(segments - 1);
    uword current_end = n - 1;
    int current_seg = segments - 1;
    while (current_seg >= 0 && current_seg < max_segments - 1) {
      int bp = last_break(current_seg, current_end);
      if (bp <= 0) break;
      bps[current_seg] = bp;
      current_end = bp - 1;
      --current_seg;
    }
    std::sort(bps.begin(), bps.end());
    breakpoint_sets[breaks - 1] = bps;
    ssr_values[breaks] = ssr_total;
    bic_values[breaks] = static_cast<double>(n) * std::log(ssr_total / static_cast<double>(n)) +
      base_constant + std::log(static_cast<double>(n)) * (order_val + 2) * segments;
  }

  int best_index = 0;
  double best_bic = bic_values[0];
  for (int i = 1; i <= max_breaks_val; ++i) {
    if (R_finite(bic_values[i]) && (i == 0 || bic_values[i] < best_bic)) {
      best_bic = bic_values[i];
      best_index = i;
    }
  }

  Rcpp::IntegerVector best_breakpoints;
  if (best_index > 0 && best_index <= max_breaks_val) {
    best_breakpoints = breakpoint_sets[best_index - 1];
  }

  Rcpp::List all_bp(max_breaks_val);
  for (int i = 0; i < max_breaks_val; ++i) all_bp[i] = breakpoint_sets[i];

  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("breakpoints") = best_breakpoints,
      Rcpp::Named("all_breakpoints") = all_bp,
      Rcpp::Named("SSR") = ssr_values,
      Rcpp::Named("BIC") = bic_values
  ));
}
