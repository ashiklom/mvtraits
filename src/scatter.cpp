#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat scatter(arma::mat mat) {
    int nr = mat.n_rows;
    int nc = mat.n_cols;
    arma::mat S = arma::zeros(nc, nc);
    arma::rowvec mrow(nc);
    for (int i = 0; i < nr; i++) {
        mrow = mat.row(i);
        S += mrow.t() * mrow;
    }
    return(S);
}
