#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

// Prototypes for exported functions

extern "C" SEXP rcpp_dummy_value();
extern "C" SEXP rcpp_fit_ar(SEXP z, SEXP p, SEXP demean);
extern "C" SEXP rcpp_bai_perron(SEXP data, SEXP order, SEXP interval, SEXP max_breaks);
extern "C" SEXP rcpp_baar(SEXP k, SEXP time, SEXP data, SEXP iterations, SEXP burn_in,
                         SEXP make_murder_p, SEXP percent, SEXP lambda, SEXP jump_p,
                         SEXP ar, SEXP progress, SEXP fit_storage);

// Simple compiled algorithm that returns the value 42.
extern "C" SEXP rcpp_dummy_value() {
  return Rcpp::wrap(42);
}

static const R_CallMethodDef CallEntries[] = {
  {"rcpp_dummy_value", (DL_FUNC)&rcpp_dummy_value, 0},
  {"rcpp_fit_ar", (DL_FUNC)&rcpp_fit_ar, 3},
  {"rcpp_bai_perron", (DL_FUNC)&rcpp_bai_perron, 4},
  {"rcpp_baar", (DL_FUNC)&rcpp_baar, 12},
  {NULL, NULL, 0}
};

extern "C" void R_init_bayesBreaks(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
