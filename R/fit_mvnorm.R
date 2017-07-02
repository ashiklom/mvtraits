#' @useDynLib mvtraits
#' @export
fit_mvnorm <- function(dat, niter = 5000, priors = list(), nchains = 3, parallel = TRUE,
                       autofit = FALSE, max_attempts = 10, threshold = 1.15, save_progress = NULL) {

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
                      setup)
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
        if (!is.null(save_progress)) {
            save_fname <- sprintf('%s.%03d', save_progress, attempt)
            saveRDS(curr_results, save_fname)
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
            rmcmc <- results2mcmclist(results_list, 'multi')
            gd <- coda::gelman.diag(rmcmc)[[1]][,1]
            exceed <- gd > threshold
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
                mu[[i]] <- colMeans(results_list[[i]][['mu']][sechalf,])
                Sigma_vec <- colMeans(results_list[[i]][['Sigma']][sechalf,])
                Sigma[[i]] <- lowerdiag2mat(Sigma_vec)
            }
        }
    }

    if (parallel) {
        parallel::stopCluster(cl)
    }

    return(results_list)
}
