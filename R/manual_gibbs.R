#' @export
fit_mvnorm <- function(dat, niter = 5000, 
                       mu0 = rep(0, ncol(dat)), Sigma_0 = diag(10, ncol(dat)),
                       v0 = ncol(dat), S0 = diag(ncol(dat)),
                       mu_init = rep(0, ncol(dat)), Sigma_init = diag(ncol(dat))) {

    # Precalculate certain quantities
    nparam <- ncol(dat)
    n <- nrow(dat)
    Sigma_0_inv <- solve(Sigma_0)

    # Setup storage
    mu_samp <- matrix(NA_real_, nrow = niter, ncol = nparam)
    Sigma_samp <- array(NA_real_, c(niter, nparam, nparam))

    has_missing <- any(is.na(dat))

    # If no values missing, pre-calculate more quantities
    if (!has_missing) {
        y <- dat
        ybar <- colMeans(y)
    }

    # Set initial conditions
    mu <- mu_init
    Sigma <- Sigma_init

    pb <- txtProgressBar(1, niter, style = 3)
    for (i in seq_len(niter)) {
        setTxtProgressBar(pb, i)
        Sigma_inv <- solve(Sigma)
        if (has_missing) {
            Sigma_chol <- chol(Sigma)
            y <- mvnorm_fill_missing(dat, mu, Sigma_chol)
            ybar <- colMeans(y)
        }
        mu <- draw_mu(ybar, n, Sigma_inv, mu0, Sigma_0_inv)
        Sigma <- draw_Sigma(y, mu, v0, S0)
        # Store outputs
        mu_samp[i,] <- mu
        Sigma_samp[i,,] <- Sigma
    }
    close(pb)
    result <- list(mu = mu_samp, Sigma = Sigma_samp)
    return(result)
}

#' @export
Rcpp::cppFunction(depends = 'RcppArmadillo',
                  code = '
    arma::mat mvnorm_fill_missing(arma::mat y, arma::vec mu, arma::mat Sigma_chol) {
        int nr = y.n_rows;
        int nc = y.n_cols;
        arma::mat Sigma_chol_inv = Sigma_chol.i();
        arma::mat x(nr, nc);
        arma::uvec j1(1);
        int nmiss;
        for (int j = 0; j < nc; j++) {
            arma::uvec jj = arma::regspace<arma::uvec>(0, j);
            arma::uvec ypres = arma::find_finite(y.col(j));
            arma::uvec ymiss = arma::find_nonfinite(y.col(j));
            nmiss = ymiss.size();
            j1.fill(j);
            x(ymiss, j1) = arma::randn(nmiss);
            y(ymiss, j1) = x(ymiss, jj) * Sigma_chol(jj, j1);
            x(ypres, j1) = y(ypres, jj) * Sigma_chol_inv(jj, j1);
        }
        return(y);
    }
                  ')

#mvnorm_fill_missing <- function(y, mu, Sigma_chol) {
    #Sigma_chol_inv <- solve(Sigma_chol)
    #nc <- ncol(y)
    #nr <- nrow(y)
    #x <- matrix(NA_real_, nr, nc)
    #for (j in seq_len(nc)) {
        #ymiss <- is.na(y[,j])
        #nmiss <- sum(ymiss)
        #x[ymiss, j] <- rnorm(nmiss)
        #x[!ymiss, j] <- y[!ymiss, 1:j, drop = FALSE] %*% Sigma_chol_inv[1:j, j]
        #y[ymiss, j] <- x[ymiss, 1:j, drop = FALSE] %*% Sigma_chol[1:j, j]
    #}
    #return(y)
#}

#' @export
Rcpp::cppFunction(depends = 'RcppArmadillo',
                  code = '
    arma::mat scatter(arma::mat mat) {
        int nr = mat.n_rows;
        int nc = mat.n_cols;
        arma::mat S = arma::zeros(nc, nc);
        arma::rowvec mrow(nc);
        for (int i = 0; i < nr; i++) {
            mrow = mat.row(i);
            S += mrow.t() * mrow;
        }
        return(S);
    }
                  ')

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
    ss <- sweep(x, 2, mu, "-")
    S_theta <- scatter(ss)
    S_n_inv <- S0 + S_theta
    S_n <- solve(S_n_inv)
    df_n <- v0 + nx
    Sigma_inv <- rWishart(1, df = df_n, Sigma = S_n)[,,1]
    Sigma <- solve(Sigma_inv)
    return(Sigma)
}

