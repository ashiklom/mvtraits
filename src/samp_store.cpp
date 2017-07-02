#include "mvtraits_headers.h"

arma::rowvec store_covmat(arma::mat X) {
    arma::vec xvec_full = arma::vectorise(arma::trimatl(X));
    arma::vec xvec = arma::nonzeros(xvec_full);
    return xvec.t();
}

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
