r_sample_mvnorm_hier <- function(
    niter, dat, groups,
    mu_global, Sigma_global, mu_group, Sigma_group,
    mu0_global, Sigma0_global_inv, mu0_group, Sigma0_group_inv,
    v0_global, S0_global, v0_group, S0_group,
    setup_bygroup
    ) {

    ugroup <- unique(groups)
    ngroup <- length(ugroup)
    n <- nrow(dat)
    m <- ncol(dat)
    mf <- m * (m + 1) / 2
    mg <- m * ngroup
    mfg <- mf * ngroup

    Sigma0_global <- solve(Sigma0_global_inv)

    Sigma_group_inv <- array(0, c(m, m, ngroup))
    Sigma0_group <- array(0, c(m, m, ngroup))

    for (g in seq_len(ngroup)) {
      Sigma_group_inv[,,g] <- solve(Sigma_group[,,g])
      Sigma0_group[,,g] <- solve(Sigma0_group_inv[,,g])
    }

    pb <- txtProgressBar()

    mu_global_samp <- matrix(NA_real_, niter, m)
    Sigma_global_samp <- matrix(NA_real_, niter, mf)
    mu_group_samp <- matrix(NA_real_, niter, mg)
    Sigma_group_samp <- matrix(NA_real_, niter, mfg)

    for (i in seq_len(niter)) {
      Sigma_global_inv <- solve(Sigma_global)

      # Within-group fit
      for (g in seq_len(ngroup)) {
        setup <- setup_bygroup[[g]]

        am <- setup[['all_missing']] + 1
        any_am <- length(am) > 0

        if (any_am) {
          amt <- setup[['missing_tail']] + 1
          o2 <- setup[['new_order']] + 1
          o21 <- setup[['revert_order']] + 1
          
          mu0_am <- mu0_group[g, am]
          Sigma0_am <- Sigma0_group[,,g][am, am, drop = FALSE]
          df0 <- 1 + v0_group[g] + m
        }

        mu <- mu_group[g,]
        Sigma <- Sigma_group[,,g]   ## CHECK SIGMA DIMENSIONS
        Sigma_inv <- solve(Sigma)
        dat_sub <- dat[groups == g,]

        y <- c_mvnorm_fill_missing(dat_sub, mu, Sigma, setup)
        ybar <- colMeans(y)

        mu_group[g,] <- c_draw_mu(ybar, nrow(y), Sigma_group_inv[,,g], mu0_group[g,], Sigma0_group_inv[,,g])
        Sigma_data <- c_draw_Sigma(y, mu_group[g,], v0_group[g], S0_group[,,g])

        if (any_am) {
          mu_group[g, am] <- random_mvnorm(1, mu0_am, Sigma0_am)
          Sigma_data2 <- Sigma_data[o2, o2]
          chol_data2 <- chol(Sigma_data2)

          Sigma_prior <- solve(rwishart(df0, S0_group[,,g]))
          Sigma_prior2 <- Sigma_prior[o2, o2]
          chol_prior2 <- chol(Sigma_prior2)

          chol_mod <- chol_data2
          chol_mod[, amt] <- chol_prior2[, amt]

          Sigma2 <- t(chol_mod) %*% chol_mod
          Sigma_final <- Sigma2[o21, o21]
        } else {
          Sigma_final <- Sigma_data
        }

        Sigma_group[,,g] <- Sigma_final
        Sigma_group_inv[,,g] <- solve(Sigma_final)

      }

      # Across-group fit
      ybar <- colMeans(mu_group)
      mu_global <- c_draw_mu(ybar, ngroup, Sigma_global_inv, mu0_global, Sigma0_global_inv)
      Sigma_global <- c_draw_Sigma(mu_group, mu_global, v0_global, S0_global)
      Sigma_global_inv <- solve(Sigma_global)

      # Store outputs
      mu_global_samp[i,] <- mu_global
      Sigma_global_samp[i,] <- store_covmat(Sigma)
      mu_group_samp[i,] <- c(mu_group)
      Sigma_group_samp[i,] <- store_covgrouparray(Sigma_group)
      setTxtProgressBar(pb, i / niter)
    }

    list(
      mu_global = mu_global_samp, Sigma_global = Sigma_global_samp,
      mu_group = mu_group_samp, Sigma_group = Sigma_group_samp
    )
}
