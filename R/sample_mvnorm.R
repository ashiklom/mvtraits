sample_mvnorm <- function(niter, dat, 
                          mu, Sigma, 
                          mu0, Sigma0_inv, 
                          v0, S0,
                          mu_samp, Sigma_samp,
                          setup) {
    # Re-arrange dimensions of Sigma_samp to fit C code
    c_Sigma_samp <- aperm(Sigma_samp, c(2, 3, 1))
    result <- c_sample_mvnorm(niter, dat, mu, Sigma,
                              mu0, Sigma0_inv, v0, S0,
                              mu_samp, c_Sigma_samp,
                              setup)
    # Return result in original dimensions
    result$Sigma <- aperm(result$Sigma, c(3, 1, 2))
    # Add names to results
    colnames(result$mu) <- names(mu)
    dimnames(result$Sigma) <- dimnames(Sigma_samp)
    return(result)
}
