solve <- function(x, ...) {
    out <- try(chol2inv(chol(x, ...)))
    if (class(out) == 'try-error') {
        print(x)
        stop('Solve failed')
    }
    return(out)
}

#' @export
draw_mu <- function(xbar, nx, Sigma_inv, mu0, Sigma_0_inv) {
    # mu | x, Sigma
    # ---------------
    # xbar -- Data matrix means
    # nx -- Number of rows in data matrix
    # Sigma_inv -- Current precision (inverse variance) matrix
    # mu0 -- Multivariate normal prior mean
    # Sigma_0_inv -- Prior precision (inverse variance) matrix
    A_n <- Sigma_0_inv + nx * Sigma_inv
    Sigma_n <- solve(A_n)
    Sigma_n_chol <- chol(Sigma_n)
    b_n <- Sigma_0_inv %*% mu0 + nx * Sigma_inv %*% xbar
    mu_n <- Sigma_n %*% b_n
    nparam <- length(mu0)
    # Random normal draw, based on Cholesky decomposition of 
    mu <- Sigma_n_chol %*% rnorm(nparam) + mu_n
    #mu <- mvtnorm::rmvnorm(1, mu_n, Sigma_n)[1,]
    return(mu)
}

#' @export
draw_Sigma <- function(x, mu, S0, v0) {
    # Sigma_star | x, mu
    # ------------------
    # x -- Data matrix
    # mu -- Mean vector
    # S0 -- Wishart prior matrix
    # v0 -- Wishart prior degrees of freedom
    nx <- nrow(x)
    nparam <- ncol(x)
    ss <- sweep(x, 2, mu, "-")
    S_theta <- scatter(ss)
    S_n_inv <- S0 + S_theta
    S_n <- solve(S_n_inv)
    df_n <- v0 + nx + nparam + 1
    Sigma_inv <- rWishart(1, df = df_n, Sigma = S_n)[,,1]
    Sigma <- solve(Sigma_inv)
    return(Sigma)
}

