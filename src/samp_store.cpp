#include "mvtraits_headers.h"

// [[Rcpp::export]]
arma::rowvec store_covmat(arma::mat X) {
    arma::vec xvec_full = arma::vectorise(arma::trimatl(X));
    arma::vec xvec = arma::nonzeros(xvec_full);
    return xvec.t();
}
