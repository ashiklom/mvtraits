#include "mvtraits_common.h"

// [[Rcpp::export]]
arma::rowvec store_covmat(arma::mat X) {
    unsigned int i;
    unsigned int j;
    unsigned int k;
    unsigned int mf = X.n_rows * (X.n_rows + 1) / 2;
    arma::vec xvec(mf);
    k = 0;
    for (j = 0; j < X.n_rows; j++) {
        for (i = j; i < X.n_rows; i++) {
            xvec(k) = X(i, j);
            k++;
        }
    }
    // for (i = 0; i < X.n_rows; i++) {
    // }
    // arma::vec xvec_full = arma::vectorise(arma::trimatl(X));
    // arma::vec xvec = arma::nonzeros(xvec_full);
    return xvec.t();
}

// [[Rcpp::export]]
arma::rowvec store_covgrouparray(arma::cube X) {
    unsigned int ngroup = X.n_slices;
    unsigned int m = X.n_rows;
    unsigned int mf = m * (m + 1) / 2;
    unsigned int mfg = mf * ngroup;
    unsigned int a; unsigned int b;
    arma::rowvec xvec(mfg, arma::fill::zeros);
    for (unsigned int i = 0; i < ngroup; i++) {
        a = i * mf;
        b = a + mf - 1;
        xvec(arma::span(a, b)) = store_covmat(X.slice(i));
    }
    return xvec;
}
