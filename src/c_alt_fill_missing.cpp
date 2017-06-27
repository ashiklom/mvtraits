#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat c_alt_fill_missing(arma::mat y, arma::vec mu, arma::mat Sigma) {
    int nr = y.n_rows;
    int nc = y.n_cols;
    arma::mat Sigma_chol = chol(Sigma);
    for (int i = 0; i < nr; i++) {
        arma::uvec m = arma::find_nonfinite(y.row(i));
        arma::uvec p = arma::find_finite(y.row(i));
        if (m.n_elem == 0) {
            continue;
        }
        if (p.n_elem == 0) {
            y.row(i) = arma::randn(nc) * Sigma_chol + mu;
            continue;
        }
        arma::vec mu_m = mu(m);
        arma::vec yp_diff = y.row(i) - mu(p);
        arma::mat Sigma_mm = Sigma[mv, mv];
        arma::mat Sigma_pp = Sigma(pv, pv);
        arma::mat Sigma_pp_inv = solve(Sigma_pp)
        arma::mat Sigma_mp = Sigma(pv, mv);
        Sigma_mpt = t(Sigma_mp);
        mu_fill = mu_m + Sigma_mpt * Sigma_pp_inv * yp_diff;
        Sigma_fill = Sigma_mm - Sigma_mpt * Sigma_pp_inv * Sigma_mp;
        y[]
        dat_filled[i, mv] = mvtnorm::rmvnorm(1, mu_fill, Sigma_fill);
    }
