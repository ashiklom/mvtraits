#' @export
random_mvnorm <- function(n, mu, Sigma) {
    if (is.matrix(mu)) {
        stopifnot(nrow(mu) %in% c(1, n))
        nc <- ncol(mu)
    } else {
        nc <- length(mu)
    }
    stopifnot(ncol(Sigma) == nc)
    Sigma_chol <- chol(Sigma)
    rand <- matrix(rnorm(n * nc), nrow = n, ncol = nc)
    out_0 <- rand %*% Sigma_chol     # (n x nc) = (n x nc) x (nc x nc)
    if (!is.matrix(mu) || nrow(mu) == 1) {
        out <- sweep(out_0, 2, mu, "+")
    } else if (all(dim(mu) == dim(out_0))) {
        out <- out_0 + mu
    } else {
        stop('Dimension mismatch')
    }
    return(out)
}

#' @export
setup_missing <- function(dat) {
    missing_vec <- apply(dat, 1, function(x) which(is.na(x)))
    missing_vec_unique <- unique(missing_vec)
    missing_pattern <- sapply(missing_vec, paste, collapse = '_')
    missing_pattern_unique <- unique(missing_pattern)
    pattern_inds <- sapply(missing_pattern_unique, function(x) which(missing_pattern == x))
    present_vec <- apply(dat, 1, function(x) which(!is.na(x)))
    present_vec_unique <- unique(present_vec)
    n_pattern <- length(missing_pattern_unique)
    out <- list(npatt = n_pattern, mlist = missing_vec_unique, plist = present_vec_unique,
                indlist = pattern_inds)
    return(out)
}

#' @export
alt_fill_missing <- function(dat, mu, Sigma, setup = NULL) {
    # Split data into missing pieces
    if (is.null(setup)) {
        setup <- setup_missing(dat)
    }
    for (i in seq_len(setup$npatt)) {
        # Rows containing that missingness pattern
        rows <- setup$indlist[[i]]
        nrows <- length(rows)
        # Columns with missing data
        m <- setup$mlist[[i]]
        if (length(m) == 0) {
            # All values present. Proceed to next pattern
            next
        }
        # Columns with present data
        p <- setup$plist[[i]]
        if (length(p) == 0) {
            # All values absent. Just do a multivariate normal draw
            dat[rows, ] <- random_mvnorm(nrows, mu, Sigma)
            next
        }
        # Partially missing. Derive parameters for partial draw.
        mu_m <- mu[m]
        mu_p <- mu[p]
        y_mup <- sweep(dat[rows, p, drop = FALSE], 2, mu_p, "-")    # (n x p)
        Sigma_mm <- Sigma[m, m, drop = FALSE]       # (m x m)
        Sigma_pp <- Sigma[p, p, drop = FALSE]       # (p x p)
        Sigma_pp_inv <- solve(Sigma_pp)             # (p x p)
        Sigma_pm <- Sigma[p, m, drop = FALSE]       # (p x m)
        Sigma_mp <- Sigma[m, p, drop = FALSE]       # (m x p)
        Sigma_prod <- Sigma_pp_inv %*% Sigma_pm     # (p x m) = (p x p) x (p x m)
        mu_0_fill <- y_mup %*% Sigma_prod           # (n x m) = (n x p) x (p x m)
        mu_fill <- sweep(mu_0_fill, 2, mu_m, "+")   # (n x m)
        Sigma_fill <- Sigma_mm - Sigma_mp %*% Sigma_prod    # (m x m) = (m x p) x (p x m)
        dat[rows, m] <- random_mvnorm(nrows, mu_fill, Sigma_fill)
    }
    return(dat)
}
