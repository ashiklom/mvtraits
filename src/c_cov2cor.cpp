#include "mvtraits_headers.h"

//[[Rcpp::export]]
arma::mat c_cov2cor(arma::mat m) {
    arma::mat d = diagmat(m);
    d = arma::inv_sympd(arma::sqrt(d));
    return d * m * d;
}
