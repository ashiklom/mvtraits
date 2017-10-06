#' @export
fit_mvnorm_hier <- function(dat, groups, niter = 5000, priors = list(), inits = list(), nchains = 3, parallel = TRUE,
                            autofit = FALSE, max_attempts = 10, keep_samples = Inf, threshold = 1.15,
                            save_progress = NULL) {

    stopifnot(is.matrix(dat), length(groups) == nrow(dat))

    chainseq <- seq_len(nchains)

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
    setup_bygroup <- lapply(ugroups, function(x) setup_missing(dat[igroups == x,]))

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

    # Draw initial conditions from priors
    mu_global <- list()
    Sigma_global <- list()
    mu_group <- list()
    Sigma_group <- list()
    for (n in chainseq) {
        mu_global[[n]] <- random_mvnorm(1, mu0_global, Sigma0_global)[1,]
        names(mu_global[[n]]) <- param_names
        Sigma_global[[n]] <- solve(rWishart(1, v0_global + nparam + 1, S0_global)[,,1])
        dimnames(Sigma_global[[n]]) <- list(param_names, param_names)

        mu_group[[n]] <- matrix(NA_real_, nrow = ngroup, ncol = nparam)
        dimnames(mu_group[[n]]) <- list(group_names, param_names)
        Sigma_group[[n]] <- array(NA_real_, c(ngroup, nparam, nparam))
        dimnames(Sigma_group[[n]]) <- list(group_names, param_names, param_names)
        for (i in seq_len(ngroup)) {
            mu_group[[n]][i,] <- random_mvnorm(1, mu0_group[i,], Sigma0_group[i,,])
            #Sigma_group[[n]][i,,] <- solve(rWishart(1, v0_group[i] + nparam + 1, S0_group[i,,])[,,1])
            Sigma_group[[n]][i,,] <- diag(1, nparam)
        }
    }
    default_inits <- list(mu_global = mu_global,
                          Sigma_global = Sigma_global,
                          mu_group = mu_group,
                          Sigma_group = Sigma_group)
    if (!is.null(inits)) {
        inits <- modifyList(default_inits, inits)
    } else {
        inits <- default_inits
    }

    samplefun <- function(n, inits) {
        sample_mvnorm_hier(niter, dat, igroups,
                           inits$mu_global[[n]], inits$Sigma_global[[n]],
                           inits$mu_group[[n]], inits$Sigma_group[[n]],
                           mu0_global, Sigma0_global,
                           mu0_group, Sigma0_group_inv,
                           v0_global, S0_global,
                           v0_group, S0_group,
                           setup_bygroup)
    }

    raw_samples <- run_until_converged(samplefun = samplefun,
                                       model_type = 'hier',
                                       inits = inits,
                                       nchains = nchains,
                                       parallel = parallel,
                                       max_attempts = max_attempts,
                                       save_progress = save_progress,
                                       threshold = threshold,
                                       keep_samples = keep_samples,
                                       autofit = autofit)

    return(raw_samples)
}
