#include "mvtraits_headers.h"

// [[Rcpp::export]]
arma::mat c_random_mvnorm(int n, arma::mat mu, arma::mat Sigma) {
    int m = mu.n_cols;
    arma::mat Sigma_chol = chol(Sigma);
    arma::mat rand(n, m, arma::fill::randn);
    arma::mat out = rand * Sigma_chol;
    if (mu.n_rows == 1) {
        out.each_row() += mu;
    } else {
        out += mu;
    }
    return out;
}
