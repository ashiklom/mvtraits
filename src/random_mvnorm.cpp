#include "mvtraits_zigg.h"

// [[Rcpp::export]]
arma::mat c_random_mvnorm(unsigned int n, arma::mat mu, arma::mat Sigma) {
    unsigned int m = mu.n_cols;
    arma::mat Sigma_chol = chol(Sigma);
    arma::mat rand(n, m);
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < m; j++) {
            rand(i, j) = zrnorm();
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
