#include <RcppArmadillo.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::vec;
using arma::uword;

// compute autocovariances up to lag p
static vec acov(const vec &y, uword p) {
  uword n = y.n_elem;
  vec centered = y - arma::mean(y);
  vec result(p + 1, arma::fill::zeros);
  for (uword lag = 0; lag <= p; ++lag) {
    double acc = 0.0;
    for (uword i = 0; i + lag < n; ++i) {
      acc += centered[i] * centered[i + lag];
    }
    result[lag] = acc / static_cast<double>(n);
  }
  return result;
}

extern "C" SEXP rcpp_fit_ar(SEXP z, SEXP p, SEXP demean) {
  vec y = Rcpp::as<vec>(z);
  y = arma::vectorise(y);
  y = y.elem(arma::find_finite(y));
  uword n = y.n_elem;
  if (n == 0) {
    Rcpp::stop("`z` must contain at least one finite observation");
  }
  int order = Rcpp::as<int>(p);
  bool do_demean = Rcpp::as<bool>(demean);
  if (order < 0) {
    Rcpp::stop("`p` must be non-negative.");
  }
  if (static_cast<uword>(order) >= n) {
    Rcpp::stop("series length must exceed the autoregressive order");
  }

  double mu_hat = do_demean ? arma::mean(y) : 0.0;
  vec centered = y - mu_hat;

  if (order == 0) {
    vec resid = centered;
    double sigma2 = arma::mean(arma::square(resid));
    if (!R_finite(sigma2) || sigma2 <= 0) sigma2 = std::numeric_limits<double>::epsilon();
    double loglik = -0.5 * static_cast<double>(n) * (std::log(2 * M_PI) + std::log(sigma2) + 1.0);
    return Rcpp::wrap(Rcpp::List::create(
        Rcpp::Named("loglikelihood") = loglik,
        Rcpp::Named("phiHat") = Rcpp::NumericVector(0),
        Rcpp::Named("sigsqHat") = sigma2,
        Rcpp::Named("muHat") = mu_hat,
        Rcpp::Named("res") = resid,
        Rcpp::Named("fits") = y - resid
    ));
  }

  vec gamma = acov(centered, order);
  arma::mat toeplitz(order, order);
  for (int i = 0; i < order; ++i) {
    for (int j = 0; j < order; ++j) {
      toeplitz(i, j) = gamma[std::abs(i - j)];
    }
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
  double sigma2 = arma::mean(arma::square(resid));
  if (!R_finite(sigma2) || sigma2 <= 0) sigma2 = std::numeric_limits<double>::epsilon();
  double loglik = -0.5 * static_cast<double>(n) * (std::log(2 * M_PI) + std::log(sigma2) + 1.0);

  return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("loglikelihood") = loglik,
      Rcpp::Named("phiHat") = phi,
      Rcpp::Named("sigsqHat") = sigma2,
      Rcpp::Named("muHat") = mu_hat,
      Rcpp::Named("res") = resid,
      Rcpp::Named("fits") = y - resid
  ));
}
