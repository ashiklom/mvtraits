r_draw_mu <- function(xbar, nx, Sigma_inv, mu0, Sigma_0_inv) {
    nSigma_inv <- nx * Sigma_inv
    A_n <- Sigma_0_inv + nSigma_inv  # (m x m)
    Sigma_n <- solve(A_n)   # (m x m)
    b_n <- mu0 %*% Sigma_0_inv + xbar %*% nSigma_inv    # (1 x m) * (m x m) = (1 x m)
    mu_n <- b_n %*% Sigma_n     # (1 x m) * (m x m) = (1 x m)
    mu <- random_mvnorm(1, mu_n, Sigma_n)
    return(mu)
}

r_draw_Sigma <- function(x, mu, v0, S0) {
    n <- nrow(x)
    m <- ncol(x)
    x_t <- t(x) - c(mu)   # Subtract from each row, efficiently; t(n x m) = (m x n)
    S_theta <- x_t %*% t(x_t)    # (m x n) * (n x m)
    S_n_inv <- S0 + S_theta
    S_n <- solve(S_n_inv)
    df_n <- v0 + n + m + 1
    Sigma_inv <- rwishart(df_n, S_n)
    Sigma <- solve(Sigma_inv)
    return(Sigma)
}

r_mvnorm_fill_missing <- function(dat, mu, Sigma, setup) {
    npatt <- setup[['npatt']]
    indlist <- setup[['indlist']]
    mlist <- setup[['mlist']]
    plist <- setup[['plist']]
    # all_missing <- setup[['all_missing']] + 1
    # dat[, all_missing] <- random_mvnorm(nrow(dat), mu[all_missing], Sigma[all_missing, all_missing, drop = FALSE])
    for (i in seq_len(npatt)) {
        rows <- indlist[[i]] + 1
        nrows <- length(rows)
        m <- mlist[[i]] + 1
        if (length(m) == 0) {
            # All values present. Proceed to next pattern.
            next
        }
        p <- plist[[i]] + 1
        if (length(p) == 0) {
            # All rows missing. Just do a multivariate normal draw
            dat[rows,] <- random_mvnorm(nrows, mu, Sigma)
            next
        }
        # Partially missing. Derive parameters for partial draw.
        y_mup <- dat[rows, p]
        y_mup <- t(t(y_mup) - mu[p])
        Sigma_pp <- Sigma[p, p]
        Sigma_mm <- Sigma[m, m]
        Sigma_pm <- Sigma[p, m]
        Sigma_mp <- Sigma[m, p]
        Sigma_prod <- solve(Sigma_pp) %*% Sigma_pm
        mu_fill <- y_mup %*% Sigma_prod
        mu_fill <- t(t(mu_fill) + mu[m])
        Sigma_fill <- Sigma_mm - Sigma_mp %*% Sigma_prod
        dat[rows, m] <- random_mvnorm(nrows, mu_fill, Sigma_fill)
    }
    return(dat)
}

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
        y <- r_mvnorm_fill_missing(dat, mu, Sigma, setup)
        ybar <- colMeans(y)

        mu <- r_draw_mu(ybar, n, Sigma_inv, mu0, Sigma0_inv)
        mu[am] <- random_mvnorm(1, mu0_am, Sigma0_am)

        Sigma_data <- r_draw_Sigma(y, mu, v0, S0)
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
