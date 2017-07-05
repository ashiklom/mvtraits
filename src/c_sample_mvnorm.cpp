#include "mvtraits_samplers.h"

//[[Rcpp::export]]
Rcpp::List c_sample_mvnorm(int niter, arma::mat dat,
                           arma::rowvec mu, arma::mat Sigma,
                           arma::rowvec mu0, arma::mat Sigma0_inv,
                           double v0, arma::mat S0,
                           Rcpp::List setup) {
    unsigned int n = dat.n_rows;
    unsigned int m = dat.n_cols;
    unsigned int mf = m * (m + 1) / 2;
    arma::mat y(arma::size(dat), arma::fill::zeros);
    arma::mat Sigma_inv = arma::inv_sympd(Sigma);
    arma::rowvec ybar(dat.n_cols, arma::fill::zeros);

    // Storage matrices
    arma::mat mu_samp(niter, m, arma::fill::zeros);
    // Store lower triangle of covariance matrix
    arma::mat Sigma_samp(niter, mf, arma::fill::zeros); 

    Progress p(niter, true);
    for (int i = 0; i < niter; i++) {
        if (i % 100 == 0) {
            if (Progress::check_abort()) {
                Rcpp::stop("Detected user interrupt.");
            }
        }
        y = c_mvnorm_fill_missing(dat, mu, Sigma, setup);
        ybar = arma::mean(y, 0);
        mu = c_draw_mu(ybar, n, Sigma_inv, mu0, Sigma0_inv);
        Sigma = c_draw_Sigma(y, mu, v0, S0);
        Sigma_inv = arma::inv_sympd(Sigma);
        mu_samp.row(i) = mu;
        Sigma_samp.row(i) = store_covmat(Sigma);
        p.increment();
    }
    return Rcpp::List::create(Rcpp::Named("mu") = mu_samp,
                              Rcpp::Named("Sigma") = Sigma_samp);
}
