#include "mvtraits_distributions.h"

//[[Rcpp::export]]
arma::rowvec c_draw_mu(arma::rowvec xbar, int nx, arma::mat Sigma_inv,
        arma::rowvec mu0, arma::mat Sigma_0_inv) {
    arma::mat nSigma_inv = nx * Sigma_inv;
    arma::mat A_n = Sigma_0_inv + nSigma_inv;       // (m x m)
    arma::mat Sigma_n = arma::inv_sympd(A_n);                 // (m x m)
    arma::mat b_n = mu0 * Sigma_0_inv + xbar * nSigma_inv;          // (1 x m) * (m x m) = (1 x m)
    arma::mat mu_n = b_n * Sigma_n;                                 // (1 x m) * (m x m) = (1 x m)
    arma::rowvec mu = c_random_mvnorm(1, mu_n, Sigma_n);
    return mu;
}
//[[Rcpp::export]]
arma::mat c_draw_Sigma(arma::mat x, arma::rowvec mu, double v0, arma::mat S0) {
    int n = x.n_rows;
    int m = x.n_cols;
    x.each_row() -= mu;             // (n x m)
    arma::mat S_theta = x.t() * x;    // (m x n) * (n x m) = (m x m)
    arma::mat S_n_inv = S0 + S_theta;
    arma::mat S_n = arma::inv_sympd(S_n_inv);
    double df_n = v0 + n + m + 1.0;
    arma::mat Sigma_inv = rwishart(df_n, S_n);
    arma::mat Sigma = arma::inv_sympd(Sigma_inv);
    return Sigma;
}
