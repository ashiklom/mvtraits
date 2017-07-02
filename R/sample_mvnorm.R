sample_mvnorm <- function(niter, dat, 
                          mu, Sigma, 
                          mu0, Sigma0_inv, 
                          v0, S0,
                          setup) {
    # Re-arrange dimensions of Sigma_samp to fit C code
    result <- c_sample_mvnorm(niter, dat, mu, Sigma,
                              mu0, Sigma0_inv, v0, S0,
                              setup)
    # Add names to results
    params <- names(mu)
    colnames(result$mu) <- params
    colnames(result$Sigma) <- lowertri_names(params)
    return(result)
}
