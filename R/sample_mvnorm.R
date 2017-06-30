sample_mvnorm <- function(niter, dat, 
                          mu, Sigma, 
                          mu0, Sigma0_inv, 
                          v0, S0,
                          mu_samp, Sigma_samp,
                          setup) {

    ndat <- nrow(dat)
    has_missing <- any(is.na(dat))
    Sigma_inv <- solve(Sigma)
    if (!has_missing) {
        y <- dat
        ybar <- colMeans(y)
    }
    pb <- txtProgressBar(1, niter, style = 3)
    for (i in seq_len(niter)) {
        setTxtProgressBar(pb, i)
        if (has_missing) {
            y <- alt_fill_missing(dat, mu, Sigma, setup)
            ybar <- colMeans(y)
        }
        mu <- draw_mu(ybar, ndat, Sigma_inv, mu0, Sigma0_inv)
        Sigma <- draw_Sigma(y, mu, v0, S0, ndat)
        if (any(abs(Sigma) > 100)) {
            stop('Sigma is too big')
        }
        Sigma_inv <- solve(Sigma)
        # Store outputs
        mu_samp[i,] <- mu
        Sigma_samp[i,,] <- Sigma
    }
    close(pb)
    result <- list(mu = mu_samp, Sigma = Sigma_samp)
    return(result)
}
