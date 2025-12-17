#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::vec;
using arma::uvec;
using arma::uword;

static double fit_segment_loglik(const vec &segment, int ar_order) {
  uword n = segment.n_elem;
  if (n == 0) return R_NegInf;
  double mean = arma::mean(segment);
  vec centered = segment - mean;
  if (ar_order == 0) {
    double sigma2 = arma::mean(arma::square(centered));
    if (!R_finite(sigma2) || sigma2 <= 0) sigma2 = std::numeric_limits<double>::epsilon();
    return -0.5 * static_cast<double>(n) * (std::log(2 * M_PI) + std::log(sigma2) + 1.0);
  }
  // Yule-Walker
  vec gamma(ar_order + 1, arma::fill::zeros);
  for (int lag = 0; lag <= ar_order; ++lag) {
    double acc = 0.0;
    for (uword i = 0; i + lag < n; ++i) acc += centered[i] * centered[i + lag];
    gamma[lag] = acc / static_cast<double>(n);
  }
  arma::mat toeplitz(ar_order, ar_order);
  for (int i = 0; i < ar_order; ++i) {
    for (int j = 0; j < ar_order; ++j) toeplitz(i, j) = gamma[std::abs(i - j)];
  }
  vec rhs = gamma.tail(ar_order);
  vec phi = arma::solve(toeplitz, rhs);
  vec resid(n, arma::fill::zeros);
  for (uword t = 0; t < n; ++t) {
    double pred = 0.0;
    for (int j = 0; j < ar_order; ++j) {
      if (t <= static_cast<uword>(j)) break;
      pred += phi[j] * centered[t - j - 1];
    }
    resid[t] = centered[t] - pred;
  }
  double sigma2 = arma::mean(arma::square(resid));
  if (!R_finite(sigma2) || sigma2 <= 0) sigma2 = std::numeric_limits<double>::epsilon();
  return -0.5 * static_cast<double>(n) * (std::log(2 * M_PI) + std::log(sigma2) + 1.0);
}

static double segmentation_loglik(const vec &data, const arma::ivec &ends, int ar_order) {
  double sum = 0.0;
  int n_ends = static_cast<int>(ends.n_elem);
  for (int i = 1; i < n_ends; ++i) {
    int start = ends[i - 1] - 1; // convert to 0-indexed
    int stop = ends[i] - 1;
    if (stop < start) return R_NegInf;
    sum += fit_segment_loglik(data.subvec(start, stop), ar_order);
  }
  return sum;
}

static arma::ivec add_breakpoint(const arma::ivec &k_ends, int ar_order, int n) {
  int constraint = (ar_order == 1) ? 2 : (2 * ar_order - 1);
  std::vector<int> full(n);
  std::iota(full.begin(), full.end(), 1);
  std::vector<int> exclude;
  exclude.reserve(k_ends.n_elem * (2 * constraint + 1));
  for (int v : k_ends) {
    int left = std::max(1, v - constraint);
    int right = std::min(n, v + constraint);
    for (int pos = left; pos <= right; ++pos) exclude.push_back(pos);
  }
  std::sort(exclude.begin(), exclude.end());
  exclude.erase(std::unique(exclude.begin(), exclude.end()), exclude.end());
  std::vector<int> diff;
  diff.reserve(n);
  std::set_difference(full.begin(), full.end(), exclude.begin(), exclude.end(), std::back_inserter(diff));
  if (diff.empty()) return k_ends;
  int idx = diff.size() == 1
  ? diff[0]
  : diff[static_cast<int>(std::floor(::unif_rand() * diff.size()))];
  arma::ivec updated(k_ends.n_elem + 1);
  for (uword i = 0; i < k_ends.n_elem; ++i) updated[i] = k_ends[i];
  updated[updated.n_elem - 1] = idx;
  updated = arma::sort(updated);
  return updated;
}

static arma::ivec remove_breakpoint(const arma::ivec &k_ends) {
  if (k_ends.n_elem <= 2) return k_ends;
  int interior = k_ends.n_elem - 2;
  int choice = (interior == 1)
    ? 1
  : static_cast<int>(std::floor(::unif_rand() * interior)) + 1;
  arma::ivec updated(k_ends.n_elem - 1);
  int ptr = 0;
  for (uword i = 0; i < k_ends.n_elem; ++i) {
    if (static_cast<int>(i) == choice) continue;
    updated[ptr++] = k_ends[i];
  }
  return updated;
}

static arma::ivec jiggle_breakpoint(const arma::ivec &k_ends, double percent, int ar_order, int n) {
  if (k_ends.n_elem <= 2) return k_ends;
  int interior = k_ends.n_elem - 2;
  int idx = static_cast<int>(std::floor(::unif_rand() * interior)) + 1;
  int bkpt = k_ends[idx];
  int wiggle = std::max(1, static_cast<int>(std::floor(n * percent)));
  int left = std::max(1, bkpt - wiggle);
  int right = std::min(n, bkpt + wiggle);
  int constraint = (ar_order == 1) ? 2 : (2 * ar_order - 1);
  left = std::max(left, k_ends[idx - 1] + constraint);
  right = std::min(right, k_ends[idx + 1] - constraint);
  if (left >= right) return k_ends;
  int proposal = left + static_cast<int>(std::floor(::unif_rand() * (right - left + 1)));
  arma::ivec updated = k_ends;
  updated[idx] = proposal;
  updated = arma::sort(updated);
  return updated;
}

extern "C" SEXP rcpp_baar(SEXP k, SEXP time, SEXP data, SEXP iterations, SEXP burn_in,
                         SEXP make_murder_p, SEXP percent, SEXP lambda, SEXP jump_p,
                         SEXP ar, SEXP progress, SEXP fit_storage) {
  Rcpp::RNGScope rng_scope;
  Rcpp::NumericVector time_vec(time);
  Rcpp::NumericVector data_vec(data);
  int iterations_val = Rcpp::as<int>(iterations);
  int burn_val = Rcpp::as<int>(burn_in);
  double make_val = Rcpp::as<double>(make_murder_p);
  double percent_val = Rcpp::as<double>(percent);
  double lambda_val = Rcpp::as<double>(lambda);
  double jump_val = Rcpp::as<double>(jump_p);
  int ar_val = Rcpp::as<int>(ar);
  bool progress_val = Rcpp::as<bool>(progress);
  bool fit_storage_val = Rcpp::as<bool>(fit_storage);

  if (time_vec.size() != data_vec.size()) Rcpp::stop("Data and time vectors must be of equal length.");
  if (iterations_val <= 0) Rcpp::stop("`iterations` must be positive.");
  int n = data_vec.size();
  if (n < 6 * ar_val) Rcpp::stop("Data insufficient for order of AR model. Try a lower order.");

  vec y = Rcpp::as<vec>(data_vec);
  Rcpp::Nullable<Rcpp::IntegerVector> k_in(k);
  arma::ivec k_vec;
  if (k_in.isNotNull()) {
    k_vec = Rcpp::as<arma::ivec>(k_in.get());
    k_vec = k_vec.elem(arma::find(k_vec > 0));
  } else {
    k_vec = arma::ivec();
  }
  arma::ivec k_ends(k_vec.n_elem + 2);
  k_ends[0] = 1;
  for (uword i = 0; i < k_vec.n_elem; ++i) k_ends[i + 1] = k_vec[i];
  k_ends[k_ends.n_elem - 1] = n;

  auto proposal_probs = [&](const arma::ivec &ends) {
    int starting_bkpts = ends.n_elem - 2;
    int nfree = std::max(1, n - 2 * ar_val - starting_bkpts);
    double total = starting_bkpts + nfree;
    double mk = make_val * (nfree / total);
    double mr = make_val * (starting_bkpts / total);
    return std::make_pair(mk, mr);
  };

  auto loglik_fn = [&](const arma::ivec &ends) { return segmentation_loglik(y, ends, ar_val); };

  auto mh_iteration = [&](arma::ivec &ends, double make_k, double murder_k) {
    double u = ::unif_rand();
    arma::ivec proposal = ends;
    double q1 = 1.0, q2 = 1.0;
    std::string type = "move";
    int constraint = (ar_val == 1) ? 5 : ar_val * 4;

    if ((arma::max(arma::diff(ends)) >= constraint && ends.n_elem < 3) ||
        (arma::max(arma::diff(ends)) >= constraint && u <= make_k)) {
      proposal = add_breakpoint(ends, ar_val, n);
      type = "add";
    } else if (u > make_k && u <= (make_k + murder_k)) {
      proposal = remove_breakpoint(ends);
      type = "sub";
    } else {
      double move_u = ::unif_rand();
      if (move_u < jump_val) {
        proposal = remove_breakpoint(add_breakpoint(ends, ar_val, n));
        type = "move";
      } else {
        proposal = jiggle_breakpoint(ends, percent_val, ar_val, n);
        type = "jiggle";
      }
    }
    double old_loglik = loglik_fn(ends);
    double new_loglik = loglik_fn(proposal);
    double delta_bic = (-2.0 * new_loglik + std::log(static_cast<double>(n)) * (proposal.n_elem - 1) * (3 + ar_val)) -
      (-2.0 * old_loglik + std::log(static_cast<double>(n)) * (ends.n_elem - 1) * (3 + ar_val));
    double ratio = (-0.5 * delta_bic) + (std::log(q1 * R::dpois(proposal.n_elem - 2, lambda_val, false)) -
                    std::log(q2 * R::dpois(ends.n_elem - 2, lambda_val, false)));
    double u_ratio = std::log(::unif_rand());
    bool accepted = R_finite(delta_bic) && R_finite(ratio) && ratio > u_ratio;
    if (accepted) ends = proposal;
    return std::make_pair(accepted, type);
  };

  // burn-in
  auto burn_probs = proposal_probs(k_ends);
  for (int i = 0; i < burn_val; ++i) {
    mh_iteration(k_ends, burn_probs.first, burn_probs.second);
  }

  int n_cols = std::max(1, static_cast<int>(k_ends.n_elem) - 2);
  Rcpp::IntegerMatrix break_mat(iterations_val, n_cols);
  std::fill(break_mat.begin(), break_mat.end(), NA_INTEGER);
  Rcpp::NumericVector bic_vals(iterations_val);
  int accept = 0;
  int add_acc = 0, sub_acc = 0, move_acc = 0, jiggle_acc = 0;
  int add_prop = 0, sub_prop = 0, move_prop = 0, jiggle_prop = 0;

  auto probs = proposal_probs(k_ends);
  for (int iter = 0; iter < iterations_val; ++iter) {
    auto result = mh_iteration(k_ends, probs.first, probs.second);
    std::string type = result.second;
    if (type == "add") add_prop++; else if (type == "sub") sub_prop++; else if (type == "move") move_prop++; else jiggle_prop++;
    if (result.first) {
      accept++;
      if (type == "add") add_acc++; else if (type == "sub") sub_acc++; else if (type == "move") move_acc++; else jiggle_acc++;
    }
    arma::ivec interior = k_ends.subvec(1, k_ends.n_elem - 2);
    int interior_len = static_cast<int>(interior.n_elem);

    if (interior_len > break_mat.ncol()) {
      Rcpp::IntegerMatrix expanded(iterations_val, interior_len);
      std::fill(expanded.begin(), expanded.end(), NA_INTEGER);

      for (int r = 0; r <= iter; ++r) {
        for (int c = 0; c < break_mat.ncol(); ++c) {
          expanded(r, c) = break_mat(r, c);
        }
      }
      break_mat = expanded;
    }

    int ncol = break_mat.ncol();
    for (int j = 0; j < ncol; ++j) {
      break_mat(iter, j) = (j < interior_len) ? interior[j] : NA_INTEGER;
    }
    double current_loglik = loglik_fn(k_ends);
    bic_vals[iter] = -2.0 * current_loglik + std::log(static_cast<double>(n)) * (k_ends.n_elem - 1) * (3 + ar_val);
  }

  Rcpp::List proposed = Rcpp::List::create(
    Rcpp::Named("add") = add_prop,
    Rcpp::Named("sub") = sub_prop,
    Rcpp::Named("move") = move_prop,
    Rcpp::Named("jiggle") = jiggle_prop
  );
  Rcpp::List accepted = Rcpp::List::create(
    Rcpp::Named("add") = add_acc,
    Rcpp::Named("sub") = sub_acc,
    Rcpp::Named("move") = move_acc,
    Rcpp::Named("jiggle") = jiggle_acc
  );

  int ncol = break_mat.ncol();
  Rcpp::IntegerVector num_bkpts(iterations_val, 0);
  for (int i = 0; i < iterations_val; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (break_mat(i, j) > 0) {
        num_bkpts[i] += 1;
      }
    }
  }

  Rcpp::List final_list = Rcpp::List::create(
    Rcpp::Named("AcceptRate") = static_cast<double>(accept) / static_cast<double>(iterations_val),
    Rcpp::Named("ProposedSteps") = proposed,
    Rcpp::Named("AcceptedSteps") = accepted,
    Rcpp::Named("BIC") = bic_vals,
    Rcpp::Named("Breakpoints") = break_mat,
    Rcpp::Named("NumBkpts") = num_bkpts
  );

  (void)progress_val; // maintained for interface parity
  (void)fit_storage_val;
  return Rcpp::wrap(final_list);
}
