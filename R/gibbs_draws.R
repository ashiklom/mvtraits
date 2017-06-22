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
    return(mu)
}

#' @export
draw_Sigma <- function(x, mu, v0, S0) {
    # Sigma_star | x, mu
    # ------------------
    # x -- Data matrix
    # mu -- Mean vector
    # v0 -- Wishart prior degrees of freedom
    # S0 -- Wishart prior scale matrix
    nx <- nrow(x)
    nparam <- ncol(x)
    xbar <- colMeans(x)
    Q0 <- cov(x) * n
    xbm0 <- (xbar - mu0) %*% t(xbar - mu0)
    Sinv_n <- S0 + Q0 + (k0 * n) / k_n * xbm0
    S_n <- solve(Sinv_n)
    df_n <- v0 + nx + nparam + 1
    Sigma_inv <- rWishart(1, df_n, S_n)[,,1]
    Sigma <- solve(Sigma_inv)
    return(Sigma)
}

