#' Check convergence of MCMC samples
#'
#' @param results_list MCMC output of `fit_mvnorm` function
#' @param model_type `character(1)`, either `"multi"` or `"hier"`
#' @param threshold Threshold for convergence of Gelman diagnostic.
#'
#' @return `logical(1)`, TRUE if converged, FALSE if not or failed to calculate
check_convergence <- function(results_list, model_type, threshold) {
    rmcmc <- results2mcmclist(results_list, model_type)
    gd <- try(coda::gelman.diag(rmcmc)[[1]][,1])
    if (class(gd) == 'try-error') {
        warning("Unable to calculate gelman diagnostic. Assuming no convergence")
        return(FALSE)
    }
    bad <- !is.finite(gd)
    if (any(bad)) {
        warning("The following parameters had non-finite convergence values: ")
        print(gd[bad])
        return(FALSE)
    }
    exceed <- gd > threshold
    converged <- all(!exceed)
    if (!converged) {
        print('The following parameters have not converged: ')
        print(gd[exceed])
    } else {
        print('All parameters have converged')
        converged <- TRUE
    }
    return(converged)
}

run_chains <- function(samplefun, inits, nchains, parallel = FALSE, cl = NULL) {
    chainseq <- seq_len(nchains)
    if (parallel) {
        parallel::clusterSetRNGStream(cl)
        curr_results <- parallel::parLapply(cl = cl, X = chainseq, fun = samplefun, inits = inits)
    } else {
        curr_results <- lapply(chainseq, samplefun, inits = inits)
    }
    return(curr_results)
}

run_until_converged <- function(samplefun,
                                model_type,
                                inits,
                                nchains,
                                parallel,
                                max_attempts,
                                save_progress,
                                threshold,
                                keep_samples,
                                autofit) {
    stopifnot(model_type %in% c('multi', 'hier'))
    if (model_type == 'multi') {
        param_names <- colnames(inits[['mu']][[1]])
    } else if (model_type == 'hier') {
        param_names <- colnames(inits[['mu_global']][[1]])
        group_names <- rownames(inits[['mu_group']][[1]])
    }
    if (parallel) {
        ncores <- min(parallel::detectCores() - 1, nchains)
        cl <- parallel::makeCluster(ncores, "FORK")
        on.exit(parallel::stopCluster(cl))
    } else {
        cl <- NULL
    }
    handle_interrupt <- function(e) {
        message("Caught user interrupt. Returning last stored results.")
        TRUE
    }
    converged <- FALSE
    attempt <- 0
    interrupted <- FALSE
    while (!converged) {
        attempt <- attempt + 1
        interrupted <- tryCatch({
            curr_results <- run_chains(samplefun = samplefun,
                                       inits = inits,
                                       nchains = nchains,
                                       parallel = parallel,
                                       cl = cl)
            FALSE   # Value for interruption
        }, interrupt = handle_interrupt)
        if (interrupted) {
            if (!exists("results_list")) {
                stop(
                    "`results_list` not found. ",
                    "Most likely, `run_chains` hasn't finished ",
                    "a full iteration yet."
                )
            }
            return(results_list)
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
            converged <- check_convergence(results_list, model_type, threshold)
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
            curr_niter <- nrow(results_list[[1]][[1]])
            start <- pmax(curr_niter - keep_samples, 1)
            keep_seq <- seq(start, curr_niter)
            prev_results <- rapply(results_list, function(x) x[keep_seq, , drop = FALSE], how = 'replace')
            for (i in seq_len(nchains)) {
                if (model_type == 'multi') {
                    # sechalf <- seq(floor(niter * 0.75), niter)
                    # mu[[i]] <- colMeans(results_list[[i]][['mu']][sechalf,])
                    # Sigma_vec <- colMeans(results_list[[i]][['Sigma']][sechalf,])
                    inits$mu[[i]] <- results_list[[i]][['mu']][curr_niter,]
                    Sigma_vec <- results_list[[i]][['Sigma']][curr_niter,]
                    inits$Sigma[[i]] <- lowerdiag2mat(Sigma_vec)
                } else if (model_type == 'hier') {
                    # sechalf <- seq(floor(niter * 0.75), niter)
                    # mu_global[[i]] <- colMeans(results_list[[i]][['mu_global']][sechalf,])
                    inits$mu_global[[i]] <- results_list[[i]][['mu_global']][curr_niter,]
                    # Sigma_global_vec <- colMeans(results_list[[i]][['Sigma_global']][sechalf,])
                    Sigma_global_vec <- results_list[[i]][['Sigma_global']][curr_niter,]
                    inits$Sigma_global[[i]] <- lowerdiag2mat(Sigma_global_vec)
                    nparam <- ncol(results_list[[i]][['mu_global']])
                    param_names <- colnames(results_list[[i]][['mu_global']])
                    ngroup <- ncol(results_list[[i]][['mu_group']]) / nparam
                    stopifnot(ngroup %% 1 == 0)
                    # mu_group_vec <- colMeans(results_list[[i]][['mu_group']][sechalf,])
                    mu_group_vec <- results_list[[i]][['mu_group']][curr_niter,]
                    inits$mu_group[[i]] <- matrix(mu_group_vec, ngroup, nparam)
                    colnames(inits$mu_group[[i]]) <- param_names
                    rownames(inits$mu_group[[i]]) <- group_names
                    # Sigma_group_vec <- colMeans(results_list[[i]][['Sigma_group']][sechalf,])
                    Sigma_group_vec <- results_list[[i]][['Sigma_group']][curr_niter,]
                    nvec <- length(Sigma_group_vec) / ngroup
                    for (j in seq_len(ngroup)) {
                        b <- j * nvec
                        a <- b - nvec + 1
                        ab <- seq(a, b)
                        inits$Sigma_group[[i]][j,,] <- lowerdiag2mat(Sigma_group_vec[ab], hier = TRUE)
                    }
                } # end model_type switch
            } # end loop over chains
        } # end assign inits
    } # end while not converged

    return(results_list)
} # end function `run_until_converged``
