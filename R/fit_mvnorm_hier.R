#' @export
fit_mvnorm_hier <- function(dat, groups, niter = 5000, priors = list()) {

    stopifnot(is.matrix(dat), length(groups) == nrow(dat))

    nparam <- ncol(dat)
    param_names <- colnames(dat)
    if (is.null(param_names)) {
        param_names <- sprintf('par%02d', seq_len(nparam))
    }

    ngroup <- length(unique(groups))
    if (is.character(groups)) {
        groups <- factor(groups)
    }
    if (is.factor(groups)) {
        group_names <- levels(groups)
    } else {
        group_names <- sprintf('group%02d', seq_len(ngroup))
    }
    igroups <- as.integer(groups)
    ugroups <- sort(unique(igroups))

    ndat <- nrow(dat)

    # Where missing, use default priors
    default_priors <- gibbs_default_priors(nparam, ngroup)
    if (!is.null(priors)) {
        priors <- modifyList(default_priors, priors)
    } else {
        priors <- default_priors
    }

    # Set priors in environment
    mu0_global <- priors[['mu_global']]
    Sigma0_global <- priors[['Sigma_global']]
    v0_global <- priors[['v_global']]
    S0_global <- priors[['S_global']]

    mu0_group <- priors[['mu_group']]
    Sigma0_group <- priors[['Sigma_group']]
    v0_group <- priors[['v_group']]
    S0_group <- priors[['S_group']]

    # Precalculate certain quantities
    Sigma0_global_inv <- solve(Sigma0_global)
    Sigma0_group_inv <- vapply(seq_len(ngroup), function(x) solve(Sigma0_group[x,,]), Sigma0_global_inv)
    Sigma0_group_inv <- aperm(Sigma0_group_inv, c(3, 1, 2))

    # Setup storage
    mu_global_samp <- matrix(NA_real_, nrow = niter, ncol = nparam)
    Sigma_global_samp <- array(NA_real_, c(niter, nparam, nparam))
    mu_group_samp <- array(NA_real_, c(niter, ngroup, nparam))
    Sigma_group_samp <- array(NA_real_, c(niter, ngroup, nparam, nparam))

    # Draw initial conditions from priors
    mu_global <- mvtnorm::rmvnorm(1, mu0_global, Sigma0_global)[1,]
    Sigma_global <- solve(rWishart(1, v0_global, S0_global)[,,1])

    mu_group <- mu_group_samp[1,,]
    Sigma_group <- Sigma_group_samp[1,,,]
    for (i in seq_len(ngroup)) {
        mu_group[i,] <- mvtnorm::rmvnorm(1, mu0_group[i,], Sigma0_group[i,,])
        Sigma_group[i,,] <- solve(rWishart(1, v0_group[i], S0_group[i,,])[,,1])
    } 

    # Set names for everything
    dimnames(mu_global_samp) <- list(NULL, param_names)
    names(mu_global) <- param_names
    dimnames(Sigma_global_samp) <- list(NULL, param_names, param_names)
    dimnames(Sigma_global) <- list(param_names, param_names)
    dimnames(mu_group_samp) <- list(NULL, group_names, param_names)
    dimnames(mu_group) <- list(group_names, param_names)
    dimnames(Sigma_group_samp) <- list(NULL, group_names, param_names, param_names)
    dimnames(Sigma_group) <- list(group_names, param_names, param_names)

    pb <- txtProgressBar(1, niter, style = 3)
    for (i in seq_len(niter)) {
        setTxtProgressBar(pb, i)
        Sigma_global_inv <- solve(Sigma_global)
        # Sample groups
        for (g in ugroups) {
            Sigma_inv <- solve(Sigma_group[g,,])
            Sigma_chol <- chol(Sigma_group[g,,])
            y <- mvnorm_fill_missing(dat[igroups == g,], mu_group[g,], Sigma_chol)
            ybar <- colMeans(y)
            ny <- nrow(y)
            mu_group[g,] <- draw_mu(ybar, ny, Sigma_inv, mu0_group[g,], Sigma0_group_inv[g,,])
            Sigma_group[g,,] <- draw_Sigma(y, mu_group[g,], v0_group[g], S0_group[g,,])
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
