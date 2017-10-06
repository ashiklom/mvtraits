r_sample_mvnorm <- function(niter, dat, mu, Sigma, mu0, Sigma0_inv,
                            v0, S0, setup) {
    n <- nrow(dat)
    m <- ncol(dat)
    mf <- m * (m + 1) / 2

    am <- setup[['all_missing']] + 1
    amt <- seq(1 + m - length(am), m)

    o1 <- seq_len(m)
    o2 <- c(o1[!o1 %in% am], am)
    o21 <- order(o2)

    Sigma_inv <- solve(Sigma)
    Sigma0 <- solve(Sigma0_inv)

    mu0_am <- mu0[am]
    Sigma0_am <- Sigma0[am, am, drop = FALSE]
    df0 <- 1 + v0 + m

    mu_samp <- matrix(NA_real_, niter, m)
    Sigma_samp <- matrix(NA_real_, niter, mf)

    pb <- txtProgressBar()

    for (i in seq_len(niter)) {
        y <- c_mvnorm_fill_missing(dat, mu, Sigma, setup)
        ybar <- colMeans(y)

        mu <- c_draw_mu(ybar, n, Sigma_inv, mu0, Sigma0_inv)
        mu[am] <- random_mvnorm(1, mu0_am, Sigma0_am)

        Sigma_data <- c_draw_Sigma(y, mu, v0, S0)
        Sigma_data2 <- Sigma_data[o2, o2]
        chol_data2 <- chol(Sigma_data2)

        Sigma_prior <- solve(rwishart(df0, S0))
        Sigma_prior2 <- Sigma_prior[o2, o2]
        chol_prior2 <- chol(Sigma_prior2)

        chol_mod <- chol_data2
        chol_mod[, amt] <- chol_prior2[, amt]

        Sigma2 <- t(chol_mod) %*% chol_mod
        Sigma <- Sigma2[o21, o21]

        Sigma_inv <- solve(Sigma)

        mu_samp[i,] <- mu
        Sigma_samp[i,] <- store_covmat(Sigma)
        setTxtProgressBar(pb, i / niter)
    }
    return(list(mu = mu_samp, Sigma = Sigma_samp))
}
