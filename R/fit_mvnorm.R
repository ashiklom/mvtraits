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

    has_missing <- any(is.na(dat))

    # If no values missing, pre-calculate more quantities
    if (!has_missing) {
        y <- dat
        ybar <- colMeans(y)
    }

    # Draw initial conditions from priors
    mu <- mvtnorm::rmvnorm(1, mu0, Sigma0)[1,]
    Sigma <- solve(rWishart(1, v0, S0)[,,1])

    # Set names on everything
    colnames(mu_samp) <- param_names
    names(mu) <- param_names
    dimnames(Sigma_samp) <- list(NULL, param_names, param_names)
    dimnames(Sigma) <- list(param_names, param_names)

    pb <- txtProgressBar(1, niter, style = 3)
    for (i in seq_len(niter)) {
        setTxtProgressBar(pb, i)
        Sigma_inv <- solve(Sigma)
        if (has_missing) {
            Sigma_chol <- chol(Sigma)
            y <- mvnorm_fill_missing(dat, mu, Sigma_chol)
            ybar <- colMeans(y)
        }
        mu <- draw_mu(ybar, ndat, Sigma_inv, mu0, Sigma0_inv)
        Sigma <- draw_Sigma(y, mu, v0, S0)
        # Store outputs
        mu_samp[i,] <- mu
        Sigma_samp[i,,] <- Sigma
    }
    close(pb)
    result <- list(mu = mu_samp, Sigma = Sigma_samp)
    return(result)
}
