// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mvnorm_fill_missing
arma::mat mvnorm_fill_missing(arma::mat y, arma::vec mu, arma::mat Sigma_chol);
RcppExport SEXP mvtraits_mvnorm_fill_missing(SEXP ySEXP, SEXP muSEXP, SEXP Sigma_cholSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_chol(Sigma_cholSEXP);
    rcpp_result_gen = Rcpp::wrap(mvnorm_fill_missing(y, mu, Sigma_chol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"mvtraits_mvnorm_fill_missing", (DL_FUNC) &mvtraits_mvnorm_fill_missing, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvtraits(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
