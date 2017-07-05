#include "mvtraits_distributions.h"

// [[Rcpp::export]]
arma::mat c_mvnorm_fill_missing (arma::mat dat, arma::rowvec mu, arma::mat Sigma, Rcpp::List setup) {
    int npatt = setup["npatt"];
    Rcpp::List indlist = setup["indlist"];
    Rcpp::List mlist = setup["mlist"];
    Rcpp::List plist = setup["plist"];
    for (int i = 0; i < npatt; i++) {
        arma::uvec rows = indlist[i];
        int nrows = rows.n_elem;
        arma::uvec m = mlist[i];
        if (m.n_elem == 0) {
            // All values present. Proceed to next pattern.
            continue;
        }
        arma::uvec p = plist[i];
        if (p.n_elem == 0) {
            // All values missing. Just do a multivariate normal draw.
            dat.rows(rows) = c_random_mvnorm(nrows, mu, Sigma);
            continue;
        }
        // Partially missing. Derive parameters for partial draw
        arma::mat y_mup = dat(rows, p);             // (rows x p)
        y_mup.each_row() -= mu.cols(p);             // (rows x p)
        arma::mat Sigma_pp = Sigma(p, p);           // (p x p)
        arma::mat Sigma_mm = Sigma(m, m);           // (m x m)
        arma::mat Sigma_pm = Sigma(p, m);           // (p x m)
        arma::mat Sigma_mp = Sigma(m, p);           // (m x p)
        arma::mat Sigma_prod = Sigma_pp.i() * Sigma_pm;     // (p x p) * (p x m) = (p x m)
        arma::mat mu_fill = y_mup * Sigma_prod;             // (rows x p) * (p x m) = (rows x m)
        mu_fill.each_row() += mu.cols(m);
        arma::mat Sigma_fill = Sigma_mm - Sigma_mp * Sigma_prod;        // (m x p) * (p x m) = (m x m)
        dat(rows, m) = c_random_mvnorm(nrows, mu_fill, Sigma_fill);
    }
    return dat;
}
