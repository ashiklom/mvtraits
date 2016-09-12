data {
    int<lower=0> Nrow;
    int<lower=0> Npar;
    int<lower=0> Ndata;
    int<lower=0> Npft;
    vector[Ndata] y;
    int indvec[Ndata];
    int sizevec[Nrow];
    int pftvec[Nrow];

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
    vector[Npar] mu_global;
    cholesky_factor_corr[Npar] L_Omega_global;
    vector<lower=0>[Npar] sigma_vec_global;

    vector[Npar] mu_pft[Npft];
    cholesky_factor_corr[Npar] L_Omega_pft[Npft];
    vector<lower=0>[Npar] sigma_vec_pft[Npft];
}

model {
    int a;
    int b;
    matrix [Npar, Npar] L_Sigma_global;
    matrix [Npar, Npar] L_Sigma_pft[Npft];

    // Priors
    mu_global ~ multi_normal(mu0, Sigma0);

    L_Omega_global ~ lkj_corr_cholesky(lkj_eta);
    sigma_vec_global ~ cauchy(cauchy_location, cauchy_scale);
    L_Sigma_global = diag_pre_multiply(sigma_vec_global, L_Omega_global);

    mu_pft ~ multi_normal_cholesky(mu_global, L_Sigma_global);

    for (i in 1:Npft) {
        L_Omega_pft[i] ~ lkj_corr_cholesky(lkj_eta);
        sigma_vec_pft[i] ~ cauchy(cauchy_location, cauchy_scale);
        L_Sigma_pft[i] = diag_pre_multiply(sigma_vec_pft[i], L_Omega_pft[i]);
    }

    /* Ragged array trick for partially missing values
     * See STAN (version 2.12.0) manual, section 13.2 (Ragged Data Structures),
     * pg. 150-151 
     */
    a = 1;
    for(i in 1:Nrow){
        b = a + sizevec[i] - 1;
        y[a:b] ~ multi_normal_cholesky(mu_pft[pftvec[i], indvec[a:b]],
                                       L_Sigma_pft[pftvec[i], indvec[a:b], indvec[a:b]]);
        a = b;
    }

}

generated quantities {
    corr_matrix[Npar] Omega_global;
    cov_matrix[Npar] Sigma_global;

    corr_matrix[Npar] Omega_pft[Npft];
    cov_matrix[Npar] Sigma_pft[Npft];

    Omega_global = multiply_lower_tri_self_transpose(L_Omega_global);
    Sigma_global = quad_form_diag(Omega_global, sigma_vec_global);

    for (i in 1:Npft) {
        Omega_pft[i] = multiply_lower_tri_self_transpose(L_Omega_pft[i]);
        Sigma_pft[i] = quad_form_diag(Omega_pft[i], sigma_vec_pft[i]);
    }
}
