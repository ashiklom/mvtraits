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
    out <- rand %*% Sigma_chol + mu
    return(out)
}

alt_fill_missing <- function(dat, mu, Sigma) {
    # Split data into missing pieces
    missing_vec <- apply(dat, 1, function(x) which(is.na(x)))
    missing_vec_unique <- unique(missing_vec)
    missing_label <- sapply(missing_vec, paste, collapse = '_')
    missing_label_unique <- unique(missing_label)
    pattern_inds <- sapply(missing_label_unique, function(x) which(missing_label == x))
    present_vec <- apply(dat, 1, function(x) which(!is.na(x)))
    present_vec_unique <- unique(present_vec)
    n <- nrow(dat)
    n_pattern <- length(missing_label_unique)
    dat_filled <- dat
    for (i in seq_len(n_pattern)) {
        # Rows containing that missingness pattern
        rows <- pattern_inds[[i]]
        nrows <- length(rows)
        # Columns with missing data
        m <- missing_vec_unique[[i]]
        if (length(m) == 0) {
            # All values present. Proceed to next pattern
            next
        }
        # Columns with present data
        p <- present_vec_unique[[i]]
        if (length(p) == 0) {
            # All values absent. Just do a multivariate normal draw
            dat_filled[rows, ] <- random_mvnorm(nrows, mu, Sigma)
            next
        }
        # Partially missing. Derive parameters for partial draw.
        mu_m <- mu[m]
        mu_p <- mu[p]
        y_mup <- t(dat[rows, p, drop = FALSE] - mu_p)
        Sigma_mm <- Sigma[m, m, drop = FALSE]
        Sigma_pp <- Sigma[p, p, drop = FALSE]
        Sigma_pp_inv <- solve(Sigma_pp)
        Sigma_pm <- Sigma[p, m, drop = FALSE]
        Sigma_pmt <- t(Sigma_pm)
        Sigma_prod <- Sigma_pmt %*% Sigma_pp_inv
        mu_fill <- t(mu_m + Sigma_prod %*% y_mup)
        Sigma_fill <- Sigma_mm - Sigma_prod %*% Sigma_pm
        dat_filled[rows, m] <- random_mvnorm(nrows, mu_fill, Sigma_fill)
    }
    #dat_filled <- dat
    #for (i in seq_len(n)) {
        #mv <- missing_vec[[i]]
        #if (length(mv) == 0) next
        #pv <- present_vec[[i]]
        #if (length(pv) == 0) {
            #dat_filled[i, ] <- mvtnorm::rmvnorm(1, mu, Sigma)
            #next
        #}
        #mu_m <- mu[mv]
        #mu_p <- mu[pv]
        #yp_diff <- t(dat[i, pv, drop = FALSE] - mu_p)
        #Sigma_mm <- Sigma[mv, mv, drop = FALSE]
        #Sigma_pp <- Sigma[pv, pv, drop = FALSE]
        #Sigma_pp_inv <- solve(Sigma_pp)
        #Sigma_mp <- Sigma[pv, mv, drop = FALSE]
        #Sigma_mpt <- t(Sigma_mp)
        #mu_fill <- mu_m + Sigma_mpt %*% Sigma_pp_inv %*% yp_diff
        #Sigma_fill <- Sigma_mm - Sigma_mpt %*% Sigma_pp_inv %*% Sigma_mp
        #dat_filled[i, mv] <- mvtnorm::rmvnorm(1, mu_fill, Sigma_fill)[1,]
    #}
    return(dat_filled)
}
