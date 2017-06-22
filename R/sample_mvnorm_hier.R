sample_mvnorm_hier <- function(niter, dat, groups,
                               mu_global, Sigma_global,
                               mu_group, Sigma_group,
                               mu0_global, Sigma0_global_inv,
                               mu0_group, Sigma0_group_inv,
                               v0_global, S0_global,
                               v0_group, S0_group,
                               mu_global_samp, Sigma_global_samp,
                               mu_group_samp, Sigma_group_samp
                               ) {

    ugroups <- sort(unique(groups))
    ngroup <- length(ugroups)

    print(Sigma_group)

    pb <- txtProgressBar(1, niter, style = 3)
    for (i in seq_len(niter)) {
        setTxtProgressBar(pb, i)
        Sigma_global_inv <- solve(Sigma_global)
        # Sample groups
        for (g in ugroups) {
            Sigma_inv <- solve(Sigma_group[g,,])
            Sigma_chol <- chol(Sigma_group[g,,])
            y <- mvnorm_fill_missing(dat[groups == g,], mu_group[g,], Sigma_chol)
            if (any(y > 1e4)) {
                print('head y')
                print(head(y))
                print('colmeans y')
                print(colMeans(y))
                print('mu')
                print(mu_group[g,])
                print('Sigma_chol')
                print(Sigma_chol)
                stop('y seems abnormally large')
            }
            ybar <- colMeans(y)
            ny <- nrow(y)
            mu_group[g,] <- draw_mu(ybar, ny, Sigma_inv, mu0_group[g,], Sigma0_group_inv[g,,])
            Sigma_group[g,,] <- draw_Sigma(y, mu_group[g,], v0_group[g], S0_group[g,,])
            if (any(Sigma_group > 1000)) {
                print('colmeans y')
                print(colMeans(y))
                print('priors')
                print(v0_group[g])
                print(S0_group[g,,])
                print('sigma')
                print(Sigma_group[g,,])
                stop('Sigma seems too big')
            }
        }
        # Sample global means
        ybar_global <- colMeans(mu_group)
        mu_global <- draw_mu(ybar_global, ngroup, Sigma_global_inv, mu0_global, Sigma0_global_inv)
        Sigma_global <- draw_Sigma(mu_group, mu_global, v0_global, S0_global)

        # Store outputs
        mu_global_samp[i,] <- mu_global
        Sigma_global_samp[i,,] <- Sigma_global
        mu_group_samp[i,,] <- mu_group
        Sigma_group_samp[i,,,] <- Sigma_group
    }
    close(pb)
    result <- list(mu_global = mu_global_samp, 
                   Sigma_global = Sigma_global_samp,
                   mu_group = mu_group_samp,
                   Sigma_group = Sigma_group_samp)
    return(result) 
}
