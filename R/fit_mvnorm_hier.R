#' @export
fit_mvnorm_hier <- function(dat, groups, niter = 5000, priors = list(), nchains = 3, parallel = TRUE,
                            autofit = FALSE, max_attempts = 10, threshold = 1.15) {

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
    setup_bygroup <- lapply(ugroups, function(x) setup_missing(dat[groups == x,]))

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
    dimnames(mu_global_samp) <- list(NULL, param_names)
    Sigma_global_samp <- array(NA_real_, c(niter, nparam, nparam))
    dimnames(Sigma_global_samp) <- list(NULL, param_names, param_names)
    mu_group_samp <- array(NA_real_, c(niter, ngroup, nparam))
    dimnames(mu_group_samp) <- list(NULL, group_names, param_names)
    Sigma_group_samp <- array(NA_real_, c(niter, ngroup, nparam, nparam))
    dimnames(Sigma_group_samp) <- list(NULL, group_names, param_names, param_names)

    # Draw initial conditions from priors
    mu_global <- list()
    Sigma_global <- list()
    mu_group <- list()
    Sigma_group <- list()
    for (n in chainseq) {
        mu_global[[n]] <- mvtnorm::rmvnorm(1, mu0_global, Sigma0_global)[1,]
        names(mu_global[[n]]) <- param_names
        Sigma_global[[n]] <- solve(rWishart(1, v0_global + nparam + 1, S0_global)[,,1])
        dimnames(Sigma_global[[n]]) <- list(param_names, param_names)

        mu_group[[n]] <- mu_group_samp[1,,]
        dimnames(mu_group[[n]]) <- list(group_names, param_names)
        Sigma_group[[n]] <- Sigma_group_samp[1,,,]
        dimnames(Sigma_group[[n]]) <- list(group_names, param_names, param_names)
        for (i in seq_len(ngroup)) {
            mu_group[[n]][i,] <- mvtnorm::rmvnorm(1, mu0_group[i,], Sigma0_group[i,,])
            #Sigma_group[[n]][i,,] <- solve(rWishart(1, v0_group[i] + nparam + 1, S0_group[i,,])[,,1])
            Sigma_group[[n]][i,,] <- diag(1, nparam)
        } 
    }

    samplefun <- function(n) {
        sample_mvnorm_hier(niter, dat, igroups,
                           mu_global[[n]], Sigma_global[[n]],
                           mu_group[[n]], Sigma_group[[n]],
                           mu0_global, Sigma0_global,
                           mu0_group, Sigma0_group_inv,
                           v0_global, S0_global,
                           v0_group, S0_group,
                           mu_global_samp, Sigma_global_samp,
                           mu_group_samp, Sigma_group_samp,
                           setup_bygroup)
    }

    if (parallel) {
        ncores <- min(parallel::detectCores() - 1, nchains)
        cl <- parallel::makeCluster(ncores, "FORK")
        parallel::clusterSetRNGStream(cl)
    }

    converged <- FALSE
    attempt <- 0
    while(!converged) {
        attempt <- attempt + 1
        if (parallel) {
            curr_results <- parallel::parLapply(cl = cl, X = chainseq, fun = samplefun)
        } else {
            curr_results <- lapply(chainseq, samplefun)
        }
        if (exists('prev_results')) {
            results_list <- combine_results(prev_results, curr_results)
        } else {
            results_list <- curr_results
        }
        if (!(nchains > 1)) {
            warning('Unable to check convergence because only one chain available.')
            converged <- TRUE
        } else {
            rmcmc <- results2mcmclist(results_list, chain2matrix_hier)
            gd <- coda::gelman.diag(rmcmc)[[1]][,1]
            exceed <- gd > 1.15
            converged <- all(!exceed)
            if (!converged) {
                print('The following parameters have not converged: ')
                print(gd[exceed])
            } else {
                print('All parameters have converged')
                converged <- TRUE
            }
        }
        if (!autofit) {
            converged <- TRUE
        }
        if (attempt >= max_attempts) {
            print(paste('Number of attempts', attempt, 
                        'exceeds max attempts', max_attempts, 
                        'but still no convergence. Returning samples as is.'))
            converged <- TRUE
        }
        if (!converged) {
            print('Resuming sampling.')
            prev_results <- results_list
            for (i in seq_len(nchains)) {
                sechalf <- seq(floor(niter * 0.75), niter)
                mu_global[[i]] <- colMeans(results_list[[i]][['mu_global']][sechalf,])
                Sigma_global[[i]] <- apply(results_list[[i]][['Sigma_global']][sechalf,,], 2:3, mean)
                mu_group[[i]] <- apply(results_list[[i]][['mu_group']][sechalf,,], 2:3, mean)
                Sigma_group[[i]] <- apply(results_list[[i]][['Sigma_group']][sechalf,,,], 2:4, mean)
            }
        }
    }

    if (parallel) {
        parallel::stopCluster(cl)
    }

    return(results_list) 
}
