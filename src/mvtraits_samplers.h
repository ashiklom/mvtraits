#ifndef _SAMPLERS_
#define _SAMPLERS_
#include "mvtraits_common.h"

arma::rowvec c_draw_mu(arma::rowvec xvar, int nx, arma::mat Sigma_inv,
        arma::rowvec mu0, arma::mat Sigma_0_inv);

arma::mat c_draw_Sigma(arma::mat x, arma::rowvec mu, double v0, arma::mat S0);

arma::mat c_mvnorm_fill_missing (arma::mat dat, arma::rowvec mu, arma::mat Sigma, Rcpp::List setup);

arma::rowvec store_covmat(arma::mat X);

arma::rowvec store_covgrouparray(arma::cube X);

#endif
