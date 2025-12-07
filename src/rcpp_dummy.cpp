#include <Rcpp.h>
#include <R_ext/Rdynload.h>

// Simple compiled algorithm that returns the value 42.

extern "C" SEXP rcpp_dummy_value();

SEXP rcpp_dummy_value() {
  return Rcpp::wrap(42);
}

static const R_CallMethodDef CallEntries[] = {
  {"rcpp_dummy_value", (DL_FUNC)&rcpp_dummy_value, 0},
  {NULL, NULL, 0}
};

extern "C" void R_init_bayesBreaks(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
