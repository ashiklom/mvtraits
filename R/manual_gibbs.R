mvnorm_fill_missing <- function(y, mu, Sigma_chol) {
    Sigma_chol_inv <- solve(Sigma_chol)
    nc <- ncol(y)
    nr <- nrow(y)
    x <- matrix(NA_real_, nr, nc)
    for (j in seq_len(nc)) {
        ymiss <- is.na(y[,j])
        nmiss <- sum(ymiss)
        x[ymiss, j] <- rnorm(nmiss)
        x[!ymiss, j] <- y[!ymiss, 1:j, drop = FALSE] %*% Sigma_chol_inv[1:j, j]
        y[ymiss, j] <- x[ymiss, 1:j, drop = FALSE] %*% Sigma_chol[1:j, j]
    }
    return(y)
}

Rcpp::cppFunction(depends = 'RcppArmadillo',
                  code = '
    arma::mat scatter(arma::mat mat) {
        int nr = mat.n_rows;
        int nc = mat.n_cols;
        arma::mat S = arma::zeros(nc, nc);
        arma::mat mrow = arma::zeros(1, nc);
        for (int i = 0; i < nr; i++) {
            mrow = mat.row(i);
            S += mrow.t() * mrow;
        }
        return(S);
    }
    
                  ')

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
    b_n <- Sigma_0_inv %*% mu0 + nx * Sigma_inv %*% xbar
    mu_n <- Sigma_n %*% b_n
    mu <- mvtnorm::rmvnorm(1, mu_n, Sigma_n)[1,]
    return(mu)
}

draw_Sigma <- function(x, mu, S0, v0) {
    # Sigma_star | x, mu
    # ------------------
    # x -- Data matrix
    # mu -- Mean vector
    # S0 -- Wishart prior matrix
    # v0 -- Wishart prior degrees of freedom
    nx <- nrow(x)
    ss <- sweep(x, 2, mu, "-")
    S_theta <- scatter(ss)
    S_n_inv <- S0 + S_theta
    S_n <- solve(S_n_inv)
    df_n <- v0 + nx
    Sigma_inv <- rWishart(1, df = df_n, Sigma = S_n)[,,1]
    Sigma <- solve(Sigma_inv)
    return(Sigma)
}

