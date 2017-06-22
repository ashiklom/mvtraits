sample_mvnorm <- function(niter, dat, 
                          mu, Sigma, 
                          mu0, Sigma0_inv, 
                          v0, S0,
                          mu_samp, Sigma_samp) {

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
        print(mu)
        print(Sigma)
        if (has_missing) {
            Sigma_chol <- chol(Sigma)
            y <- mvnorm_fill_missing(dat, mu, Sigma_chol)
            ybar <- colMeans(y)
        }
        print(ybar)
        .prompt <- readline(prompt = 'Continue: [Enter]')
        mu <- draw_mu(ybar, ndat, Sigma_inv, mu0, Sigma0_inv)
        Sigma <- draw_Sigma(y, mu, v0, S0)
        Sigma_inv <- solve(Sigma)
        # Store outputs
        mu_samp[i,] <- mu
        Sigma_samp[i,,] <- Sigma
    }
    close(pb)
    result <- list(mu = mu_samp, Sigma = Sigma_samp)
    return(result)
}
