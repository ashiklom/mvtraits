data {
    int<lower=0> Nrow;
    int<lower=0> Npar;
    int<lower=0> Ndata;
    vector[Npar] y[Nrow];
    int pres[Ndata, 2];

    // Priors
    // Mean vector (mu)
    vector[Npar] mu0;
    vector[Npar] sigma0;

    // Standard deviation vector (sigma)
    real cauchy_location;
    real cauchy_scale;
}

parameters {
    vector[Npar] mu;
    vector<lower=0>[Npar] sigma;
}

model {
    // Priors
    mu ~ normal(mu0, sigma0);
    sigma ~ cauchy(cauchy_location, cauchy_scale);
    
    // Likelihood
    for (i in 1:Ndata) {
        y[pres[i,1], pres[i,2]] ~ normal(mu0[pres[i,2]], sigma0[pres[i,2]]);
    }
}

generated quantities {
    vector[Npar] sigma2;
    sigma2 = sigma .* sigma;
}
