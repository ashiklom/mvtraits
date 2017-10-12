// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// c_cov2cor
arma::mat c_cov2cor(arma::mat m);
RcppExport SEXP _mvtraits_c_cov2cor(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(c_cov2cor(m));
    return rcpp_result_gen;
END_RCPP
}
// c_draw_mu
arma::rowvec c_draw_mu(arma::rowvec xbar, unsigned int nx, arma::mat Sigma_inv, arma::rowvec mu0, arma::mat Sigma_0_inv);
RcppExport SEXP _mvtraits_c_draw_mu(SEXP xbarSEXP, SEXP nxSEXP, SEXP Sigma_invSEXP, SEXP mu0SEXP, SEXP Sigma_0_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::rowvec >::type xbar(xbarSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_inv(Sigma_invSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_0_inv(Sigma_0_invSEXP);
    rcpp_result_gen = Rcpp::wrap(c_draw_mu(xbar, nx, Sigma_inv, mu0, Sigma_0_inv));
    return rcpp_result_gen;
END_RCPP
}
// c_draw_Sigma
arma::mat c_draw_Sigma(arma::mat x, arma::rowvec mu, double v0, arma::mat S0);
RcppExport SEXP _mvtraits_c_draw_Sigma(SEXP xSEXP, SEXP muSEXP, SEXP v0SEXP, SEXP S0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S0(S0SEXP);
    rcpp_result_gen = Rcpp::wrap(c_draw_Sigma(x, mu, v0, S0));
    return rcpp_result_gen;
END_RCPP
}
// c_sample_mvnorm_hier
Rcpp::List c_sample_mvnorm_hier(unsigned int niter, arma::mat dat, arma::uvec groups, arma::rowvec mu_global, arma::mat Sigma_global, arma::mat mu_group, arma::cube Sigma_group, arma::rowvec mu0_global, arma::mat Sigma0_global_inv, arma::mat mu0_group, arma::cube Sigma0_group_inv, double v0_global, arma::mat S0_global, arma::vec v0_group, arma::cube S0_group, Rcpp::List setup_bygroup);
RcppExport SEXP _mvtraits_c_sample_mvnorm_hier(SEXP niterSEXP, SEXP datSEXP, SEXP groupsSEXP, SEXP mu_globalSEXP, SEXP Sigma_globalSEXP, SEXP mu_groupSEXP, SEXP Sigma_groupSEXP, SEXP mu0_globalSEXP, SEXP Sigma0_global_invSEXP, SEXP mu0_groupSEXP, SEXP Sigma0_group_invSEXP, SEXP v0_globalSEXP, SEXP S0_globalSEXP, SEXP v0_groupSEXP, SEXP S0_groupSEXP, SEXP setup_bygroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu_global(mu_globalSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_global(Sigma_globalSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu_group(mu_groupSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma_group(Sigma_groupSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu0_global(mu0_globalSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma0_global_inv(Sigma0_global_invSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu0_group(mu0_groupSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type Sigma0_group_inv(Sigma0_group_invSEXP);
    Rcpp::traits::input_parameter< double >::type v0_global(v0_globalSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S0_global(S0_globalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type v0_group(v0_groupSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type S0_group(S0_groupSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type setup_bygroup(setup_bygroupSEXP);
    rcpp_result_gen = Rcpp::wrap(c_sample_mvnorm_hier(niter, dat, groups, mu_global, Sigma_global, mu_group, Sigma_group, mu0_global, Sigma0_global_inv, mu0_group, Sigma0_group_inv, v0_global, S0_global, v0_group, S0_group, setup_bygroup));
    return rcpp_result_gen;
END_RCPP
}
// c_sample_mvnorm
Rcpp::List c_sample_mvnorm(unsigned int niter, arma::mat dat, arma::rowvec mu, arma::mat Sigma, arma::rowvec mu0, arma::mat Sigma0_inv, double v0, arma::mat S0, Rcpp::List setup);
RcppExport SEXP _mvtraits_c_sample_mvnorm(SEXP niterSEXP, SEXP datSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP mu0SEXP, SEXP Sigma0_invSEXP, SEXP v0SEXP, SEXP S0SEXP, SEXP setupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma0_inv(Sigma0_invSEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type setup(setupSEXP);
    rcpp_result_gen = Rcpp::wrap(c_sample_mvnorm(niter, dat, mu, Sigma, mu0, Sigma0_inv, v0, S0, setup));
    return rcpp_result_gen;
END_RCPP
}
// c_mvnorm_fill_missing
arma::mat c_mvnorm_fill_missing(arma::mat dat, arma::rowvec mu, arma::mat Sigma, Rcpp::List setup);
RcppExport SEXP _mvtraits_c_mvnorm_fill_missing(SEXP datSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP setupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type setup(setupSEXP);
    rcpp_result_gen = Rcpp::wrap(c_mvnorm_fill_missing(dat, mu, Sigma, setup));
    return rcpp_result_gen;
END_RCPP
}
// c_random_mvnorm
arma::mat c_random_mvnorm(unsigned int n, arma::mat mu, arma::mat Sigma);
RcppExport SEXP _mvtraits_c_random_mvnorm(SEXP nSEXP, SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(c_random_mvnorm(n, mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// rwishart
arma::mat rwishart(double df, arma::mat S);
RcppExport SEXP _mvtraits_rwishart(SEXP dfSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(rwishart(df, S));
    return rcpp_result_gen;
END_RCPP
}
// store_covmat
arma::rowvec store_covmat(arma::mat X);
RcppExport SEXP _mvtraits_store_covmat(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(store_covmat(X));
    return rcpp_result_gen;
END_RCPP
}
// store_covgrouparray
arma::rowvec store_covgrouparray(arma::cube X);
RcppExport SEXP _mvtraits_store_covgrouparray(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(store_covgrouparray(X));
    return rcpp_result_gen;
END_RCPP
}
// set_R_seed
void set_R_seed();
RcppExport SEXP _mvtraits_set_R_seed() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    set_R_seed();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mvtraits_c_cov2cor", (DL_FUNC) &_mvtraits_c_cov2cor, 1},
    {"_mvtraits_c_draw_mu", (DL_FUNC) &_mvtraits_c_draw_mu, 5},
    {"_mvtraits_c_draw_Sigma", (DL_FUNC) &_mvtraits_c_draw_Sigma, 4},
    {"_mvtraits_c_sample_mvnorm_hier", (DL_FUNC) &_mvtraits_c_sample_mvnorm_hier, 16},
    {"_mvtraits_c_sample_mvnorm", (DL_FUNC) &_mvtraits_c_sample_mvnorm, 9},
    {"_mvtraits_c_mvnorm_fill_missing", (DL_FUNC) &_mvtraits_c_mvnorm_fill_missing, 4},
    {"_mvtraits_c_random_mvnorm", (DL_FUNC) &_mvtraits_c_random_mvnorm, 3},
    {"_mvtraits_rwishart", (DL_FUNC) &_mvtraits_rwishart, 2},
    {"_mvtraits_store_covmat", (DL_FUNC) &_mvtraits_store_covmat, 1},
    {"_mvtraits_store_covgrouparray", (DL_FUNC) &_mvtraits_store_covgrouparray, 1},
    {"_mvtraits_set_R_seed", (DL_FUNC) &_mvtraits_set_R_seed, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvtraits(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
