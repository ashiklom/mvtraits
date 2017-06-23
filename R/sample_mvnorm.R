sample_mvnorm <- function(niter, dat, 
                          mu, Sigma, 
                          mu0, Sigma0_inv, 
                          v0, S0,
                          mu_samp, Sigma_samp) {

    ndat <- nrow(dat)
    nmiss <- rowSums(is.na(dat))
    neff <- sum(1/nmiss^sqrt(2))
    #neff <- mean(colSums(!is.na(dat)))
    has_missing <- any(is.na(dat))
    Sigma_inv <- solve(Sigma)
    if (!has_missing) {
        y <- dat
        ybar <- colMeans(y)
    }
    pb <- txtProgressBar(1, niter, style = 3)
    #ybar <- colMeans(dat, na.rm = TRUE)
    for (i in seq_len(niter)) {
        setTxtProgressBar(pb, i)
        if (has_missing) {
            Sigma_chol <- chol(Sigma)
            y <- mvnorm_fill_missing(dat, mu, Sigma_chol)
            ybar <- colMeans(y)
            print(ybar)
        }
        #print(mu)
        #print(Sigma)
        #print(ybar)
        #.prompt <- readline(prompt = 'Continue: [Enter]')
        mu <- draw_mu(ybar, neff, Sigma_inv, mu0, Sigma0_inv)
        print(mu)
        Sigma <- draw_Sigma(y, mu, v0, S0, neff)
        print(Sigma)
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
