#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat mvnorm_fill_missing(arma::mat y, arma::vec mu, arma::mat Sigma_chol) {
    int nr = y.n_rows;
    int nc = y.n_cols;
    arma::mat Sigma_chol_inv = Sigma_chol.i();
    arma::mat x(nr, nc);
    arma::uvec j1(1);
    int nmiss;
    for (int j = 0; j < nc; j++) {
        arma::uvec jj = arma::regspace<arma::uvec>(0, j);
        arma::uvec ypres = arma::find_finite(y.col(j));
        arma::uvec ymiss = arma::find_nonfinite(y.col(j));
        nmiss = ymiss.size();
        j1.fill(j);
        x(ymiss, j1) = arma::randn(nmiss);
        y(ymiss, j1) = x(ymiss, jj) * Sigma_chol(jj, j1);
        x(ypres, j1) = y(ypres, jj) * Sigma_chol_inv(jj, j1);
    }
    return(y);
}

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
