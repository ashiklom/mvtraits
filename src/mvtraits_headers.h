#ifndef _mvtraits_headers_
#define _mvtraits_headers_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat c_random_mvnorm(int n, arma::mat mu, arma::mat Sigma);

arma::mat rwishart(double df, arma::mat S);

arma::rowvec c_draw_mu(arma::rowvec xvar, int nx, arma::mat Sigma_inv,
        arma::rowvec mu0, arma::mat Sigma_0_inv);

arma::mat c_draw_Sigma(arma::mat x, arma::rowvec mu, double v0, arma::mat S0);

arma::mat c_alt_fill_missing (arma::mat dat, arma::rowvec mu, arma::mat Sigma, Rcpp::List setup);

#endif
