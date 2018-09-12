#include "mvtraits_samplers.h"

//[[Rcpp::depends(RcppProgress)]]

//[[Rcpp::export]]
Rcpp::List c_sample_mvnorm_hier(unsigned int niter, arma::mat dat, arma::uvec groups,
                                arma::rowvec mu_global, arma::mat Sigma_global,
                                arma::mat mu_group, arma::cube Sigma_group,
                                arma::rowvec mu0_global, arma::mat Sigma0_global_inv,
                                arma::mat mu0_group, arma::cube Sigma0_group_inv,
                                double v0_global, arma::mat S0_global,
                                arma::vec v0_group, arma::cube S0_group,
                                Rcpp::List setup_bygroup) {
    arma::uvec ugroups = arma::unique(groups);
    unsigned int ngroup = ugroups.n_elem;
    unsigned int n = dat.n_rows;
    unsigned int m = dat.n_cols;
    unsigned int mf = m * (m + 1) / 2;
    unsigned int mg = m * ngroup;
    unsigned int mfg = mf * ngroup;

    unsigned int i; unsigned int g;

    // Initialize values
    arma::mat Sigma_global_inv(n, m, arma::fill::zeros);
    arma::cube Sigma_group_inv(n, m, ngroup, arma::fill::zeros);
    arma::rowvec mu(m, arma::fill::zeros);
    arma::mat Sigma(m, m, arma::fill::zeros);
    arma::mat Sigma_inv(m, m, arma::fill::zeros);
    arma::rowvec ybar(m, arma::fill::zeros);
    unsigned int ny;

    // Initialize storage
    arma::mat mu_global_samp(niter, m, arma::fill::zeros);
    arma::mat Sigma_global_samp(niter, mf, arma::fill::zeros);
    arma::mat mu_group_samp(niter, mg, arma::fill::zeros);
    arma::mat Sigma_group_samp(niter, mfg, arma::fill::zeros);

    Progress p(niter, true);

    // Sampling
    for (i = 0; i < niter; i++) {
        if (i % 100 == 0) {
            if (Progress::check_abort()) {
                Rcpp::stop("Detected user interrupt.");
            }
        }
        Sigma_global_inv = arma::inv_sympd(Sigma_global);
        // Within-group mean
        for (g = 0; g < ngroup; g++) {
            Rcpp::List setup = setup_bygroup[g];
            mu = mu_group.row(g);
            Sigma = Sigma_group.slice(g);
            Sigma_inv = arma::inv_sympd(Sigma);
            arma::mat dat_sub = dat.rows(arma::find(groups == (g + 1)));
            arma::mat y = c_mvnorm_fill_missing(dat_sub, mu, Sigma, setup);
            ybar = arma::mean(y, 0);
            mu_group.row(g) = c_draw_mu(ybar, y.n_rows, Sigma_inv, mu0_group.row(g), Sigma0_group_inv.slice(g));
            Sigma_group.slice(g) = c_draw_Sigma(y, mu_group.row(g), v0_group(g), S0_group.slice(g));
        }
        // Across group mean
        ybar = arma::mean(mu_group, 0);
        mu_global = c_draw_mu(ybar, ngroup, Sigma_global_inv, mu0_global, Sigma0_global_inv);
        Sigma_global = c_draw_Sigma(mu_group, mu_global, v0_global, S0_global);
        Sigma_global_inv = arma::inv_sympd(Sigma_global);

        // Store outputs
        mu_global_samp.row(i) = mu_global;
        Sigma_global_samp.row(i) = store_covmat(Sigma_global);
        mu_group_samp.row(i) = arma::trans(arma::vectorise(mu_group));
        Sigma_group_samp.row(i) = store_covgrouparray(Sigma_group);
        p.increment();
    }
    return Rcpp::List::create(Rcpp::Named("mu_global") = mu_global_samp,
                              Rcpp::Named("Sigma_global") = Sigma_global_samp,
                              Rcpp::Named("mu_group") = mu_group_samp,
                              Rcpp::Named("Sigma_group") = Sigma_group_samp);
}
