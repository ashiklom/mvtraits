#' @useDynLib mvtraits
#' @export
fit_mvnorm <- function(dat, niter = 5000, priors = list(), nchains = 3, parallel = TRUE,
                       repeat_until_converged = FALSE) {

    chainseq <- seq_len(nchains)

    nparam <- ncol(dat)
    param_names <- colnames(dat)
    if (is.null(param_names)) {
        param_names <- sprintf('par%02d', seq_len(nparam))
    }

    ndat <- nrow(dat)
    setup <- setup_missing(dat)

    # Where missing, use default priors
    default_priors <- gibbs_default_priors(nparam)
    if (!is.null(priors)) {
        priors <- modifyList(default_priors, priors)
    } else {
        priors <- default_priors
    }

    # Set priors in environment
    mu0 <- priors[['mu_global']]
    Sigma0 <- priors[['Sigma_global']]
    v0 <- priors[['v_global']]
    S0 <- priors[['S_global']]

    Sigma0_inv <- solve(Sigma0)

    # Setup storage
    mu_samp <- matrix(NA_real_, nrow = niter, ncol = nparam)
    colnames(mu_samp) <- param_names
    Sigma_samp <- array(NA_real_, c(niter, nparam, nparam))
    dimnames(Sigma_samp) <- list(NULL, param_names, param_names)

    # Draw initial conditions from priors
    mu <- list()
    Sigma <- list()
    for (n in chainseq) {
        mu[[n]] <- mvtnorm::rmvnorm(1, mu0, Sigma0)[1,]
        names(mu[[n]]) <- param_names
        Sigma[[n]] <- solve(rWishart(1, v0 + nparam + 1, S0)[,,1])
        dimnames(Sigma[[n]]) <- list(param_names, param_names)
    }

    samplefun <- function(n) {
        sample_mvnorm(niter, dat, mu[[n]], Sigma[[n]],
                      mu0, Sigma0_inv, v0, S0,
                      mu_samp, Sigma_samp, setup)
    }

    if (parallel) {
        ncores <- min(parallel::detectCores() - 1, nchains)
        cl <- parallel::makeCluster(ncores, "FORK")
        parallel::clusterSetRNGStream(cl)
    }

    converged <- FALSE
    while(!converged) {
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
            rmcmc <- results2mcmclist(results_list, chain2matrix_multi)
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
        if (!repeat_until_converged) {
            converged <- TRUE
        }
        if (!converged) {
            print('Resuming sampling.')
            prev_results <- results_list
            for (i in seq_len(nchains)) {
                sechalf <- seq(floor(niter * 0.75), niter)
                mu[[i]] <- colMeans(results_list[[i]][['mu']][sechalf,])
                Sigma[[i]] <- apply(results_list[[i]][['Sigma']][sechalf,,], 2:3, mean)
            }
        }
    }

    if (parallel) {
        parallel::stopCluster(cl)
    }

    return(results_list)
}
