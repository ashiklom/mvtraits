sample_mvnorm_hier <- function(niter, dat, groups,
                               mu_global, Sigma_global,
                               mu_group, Sigma_group,
                               mu0_global, Sigma0_global_inv,
                               mu0_group, Sigma0_group_inv,
                               v0_global, S0_global,
                               v0_group, S0_group,
                               mu_global_samp, Sigma_global_samp,
                               mu_group_samp, Sigma_group_samp,
                               setup_bygroup
                               ) {

    # Convert dimensions
    c_Sigma_group <- aperm(Sigma_group, c(2, 3, 1))
    c_Sigma0_group_inv <- aperm(Sigma0_group_inv, c(2, 3, 1))
    c_S0_group <- aperm(S0_group, c(2, 3, 1))
    c_Sigma_global_samp <- aperm(Sigma_global_samp, c(2, 3, 1))
    c_mu_group_samp <- aperm(mu_group_samp, c(2, 3, 1))
    c_Sigma_group_samp <- aperm(Sigma_group_samp, c(1, 3, 4, 2))
    c_Sigma_group_samp <- array2field(Sigma_group_samp)

    # Run sampler
    result <- c_sample_mvnorm_hier(niter, dat, groups,
                                   mu_global, Sigma_global,
                                   mu_group, c_Sigma_group,
                                   mu0_global, Sigma0_global_inv,
                                   mu0_group, c_Sigma0_group_inv,
                                   v0_global, S0_global,
                                   v0_group, c_S0_group,
                                   mu_global_samp, c_Sigma_global_samp,
                                   c_mu_group_samp, c_Sigma_group_samp,
                                   setup_bygroup)

    # Convert back dimensions
    result$mu_group <- aperm(result$mu_group, c(3, 1, 2))
    result$Sigma_global <- aperm(result$Sigma_global, c(3, 1, 2))
    result$Sigma_group <- aperm(field2array(result$Sigma_group), c(1, 4, 2, 3))

    # Add names
    dimnames(result$mu_global) <- dimnames(mu_global_samp)
    dimnames(result$Sigma_global) <- dimnames(Sigma_global_samp)
    dimnames(result$mu_group) <- dimnames(mu_group_samp)
    dimnames(result$Sigma_group) <- dimnames(Sigma_group_samp)
    return(result)
}

array2field <- function(a) {
    dims <- dim(a)
    ndims <- length(dims)
    stopifnot(ndims <= 4)
    cube_list <- lapply(seq_len(dims[1]), function(x) a[x,,,])
    out <- matrix(cube_list, dims[1], 1)
    return(out)
}

field2array <- function(f) {
    d1 <- nrow(f)
    d2p <- dim(f[[1]])
    dims <- c(d1, d2p)
    a <- array(NA_real_, dims)
    for (i in seq_len(d1)) {
        a[i,,,] <- f[[i]]
    }
    return(a)
}

