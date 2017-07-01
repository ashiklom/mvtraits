#include "mvtraits_headers.h"

//[[Rcpp::export]]
Rcpp::List c_sample_mvnorm(int niter, arma::mat dat,
                           arma::rowvec mu, arma::mat Sigma,
                           arma::rowvec mu0, arma::mat Sigma0_inv,
                           double v0, arma::mat S0,
                           arma::mat mu_samp, arma::cube Sigma_samp,
                           Rcpp::List setup) {
    int n = dat.n_rows;
    arma::mat y(arma::size(dat), arma::fill::zeros);
    arma::mat Sigma_inv = arma::inv_sympd(Sigma);
    arma::rowvec ybar(dat.n_cols, arma::fill::zeros);
    for (int i = 0; i < niter; i++) {
        y = c_alt_fill_missing(dat, mu, Sigma, setup);
        ybar = arma::mean(y, 0);
        mu = c_draw_mu(ybar, n, Sigma_inv, mu0, Sigma0_inv);
        Sigma = c_draw_Sigma(y, mu, v0, S0);
        Sigma_inv = arma::inv_sympd(Sigma);
        mu_samp.row(i) = mu;
        Sigma_samp.slice(i) = Sigma;
    }
    return Rcpp::List::create(Rcpp::Named("mu") = mu_samp,
                              Rcpp::Named("Sigma") = Sigma_samp);
}
