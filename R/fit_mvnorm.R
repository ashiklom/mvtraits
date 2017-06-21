#' @useDynLib mvtraits
#' @export
fit_mvnorm <- function(dat, niter = 5000, priors = list(), nchains = 3, parallel = TRUE) {

    chainseq <- seq_len(nchains)

    nparam <- ncol(dat)
    param_names <- colnames(dat)
    if (is.null(param_names)) {
        param_names <- sprintf('par%02d', seq_len(nparam))
    }

    ndat <- nrow(dat)

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
        Sigma[[n]] <- solve(rWishart(1, v0, S0)[,,1])
        dimnames(Sigma[[n]]) <- list(param_names, param_names)
    }

    ncores <- min(parallel::detectCores() - 1, nchains)
    cl <- parallel::makeCluster(ncores, "FORK")

    samplefun <- function(n) {
        sample_mvnorm(niter, dat, mu[[n]], Sigma[[n]],
                      mu0, Sigma0, v0, S0,
                      mu_samp, Sigma_samp)
    }

    results_list <- parallel::parLapply(cl = cl, X = chainseq, fun = samplefun)

    return(results_list)
}
