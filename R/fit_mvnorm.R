#' @useDynLib mvtraits
#' @export
fit_mvnorm <- function(dat, niter = 5000, priors = list()) {

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
    Sigma_samp <- array(NA_real_, c(niter, nparam, nparam))

    # If no values missing, pre-calculate more quantities
    # Draw initial conditions from priors
    mu <- mvtnorm::rmvnorm(1, mu0, Sigma0)[1,]
    Sigma <- solve(rWishart(1, v0, S0)[,,1])

    # Set names on everything
    colnames(mu_samp) <- param_names
    names(mu) <- param_names
    dimnames(Sigma_samp) <- list(NULL, param_names, param_names)
    dimnames(Sigma) <- list(param_names, param_names)

    result <- sample_mvnorm(niter, dat, mu, Sigma,
                            mu0, Sigma0, v0, S0,
                            mu_samp, Sigma_samp)
    return(result)
}
