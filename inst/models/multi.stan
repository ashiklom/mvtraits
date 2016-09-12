data {
    int<lower=0> Nrow;
    int<lower=0> Npar;
    int<lower=0> Ndata;
    vector[Ndata] y;
    int indvec[Ndata];
    int sizevec[Nrow];

    // Priors
    // Mean vector (mu)
    vector[Npar] mu0;
    cov_matrix[Npar] Sigma0;

    // Covariance matrix (Sigma)
    real cauchy_location;
    real cauchy_scale;
    real lkj_eta;
}

parameters {
    vector[Npar] mu;
    cholesky_factor_corr[Npar] L_Omega;
    vector<lower=0>[Npar] sigma_vec;
}

model {
    int a;
    int b;
    matrix [Npar, Npar] L_Sigma;

    // Priors
    mu ~ multi_normal(mu0, Sigma0);

    L_Omega ~ lkj_corr_cholesky(lkj_eta);
    sigma_vec ~ cauchy(cauchy_location, cauchy_scale);
    L_Sigma = diag_pre_multiply(sigma_vec, L_Omega);

    /* Ragged array trick for partially missing values
     * See STAN (version 2.12.0) manual, section 13.2 (Ragged Data Structures),
     * pg. 150-151 
     */
    a = 1;
    for(i in 1:Nrow){
        b = a + sizevec[i] - 1;
        y[a:b] ~ multi_normal_cholesky(mu[indvec[a:b]],
                                       L_Sigma[indvec[a:b], indvec[a:b]]);
        a = b;
    }

}

generated quantities {
    corr_matrix[Npar] Omega;
    cov_matrix[Npar] Sigma;
    Omega = multiply_lower_tri_self_transpose(L_Omega);
    Sigma = quad_form_diag(Omega, sigma_vec);
}
