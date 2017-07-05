#include "mvtraits_zigg.h"

// [[Rcpp::export]]
arma::mat c_random_mvnorm(int n, arma::mat mu, arma::mat Sigma) {
    int m = mu.n_cols;
    arma::mat Sigma_chol = chol(Sigma);
    arma::mat rand(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            rand(i, j) = zigg.norm();
        }
    }
    arma::mat out = rand * Sigma_chol;
    if (mu.n_rows == 1) {
        out.each_row() += mu;
    } else {
        out += mu;
    }
    return out;
}
