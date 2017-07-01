#ifndef _mvtraits_headers_
#define _mvtraits_headers_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat c_random_mvnorm(int n, arma::mat mu, arma::mat Sigma);

#endif
