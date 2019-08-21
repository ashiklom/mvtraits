#' @export
add_slopes <- function(results_list, hier = FALSE, ngroups = NULL) {
    vars <- names(results_list[[1]])
    for (i in seq_along(results_list)) {
        if (hier) {
            results_list[[i]][['RMAslope_global']] <- calc_slopes(results_list[[i]][['Sigma_global']])
            results_list[[i]][['RMAslope_group']] <- calc_slopes(results_list[[i]][['Sigma_group']], hier = TRUE, ngroups = ngroups)
        } else {
            results_list[[i]][['RMAslope']] <- calc_slopes(results_list[[i]][['Sigma']])
        }
    }
    return(results_list)
}

#' @export
calc_slopes <- function(sigma_samples, hier = FALSE, ngroups = NULL) {
    dims <- dim(sigma_samples)
    cdims <- dims
    if (hier && !isTRUE(is.finite(ngroups))) stop('If hierarchical, ngroups must be provided')
    ng <- ifelse(hier, ngroups, 1)
    nvec <- dims[2] / ng
    nc <- 0.5 * (sqrt(8 * nvec + 1) - 1)
    nc2 <- nc - 1
    nm <- nc2 * (nc2 + 1) / 2
    cdims[2] <- nm * ng
    cormat_samples <- array(0, cdims)
    # Determine colnames for correlation matrix
    sigma_cnames <- colnames(sigma_samples)
    names_split <- do.call(rbind, strsplit(sigma_cnames, split = varsep_esc))
    add <- ifelse(hier, 1, 0)
    ii <- 1 + add
    jj <- 2 + add
    param_names <- apply(names_split[seq_len(nc), -jj, drop = FALSE], 1, paste, collapse = varsep) 
    corr_names_mat <- names_split[names_split[,ii] != names_split[,jj], , drop = FALSE]
    corr_names <- apply(corr_names_mat, 1, paste, collapse = varsep)
    for (i in seq_len(dims[1])) {
        if (hier) {
            # Hierarchical model -- ngroup (slow) x nparam (fast)
            for (j in seq_len(ngroups)) {
                b <- j * nvec
                a <- b - nvec + 1
                ab <- seq(a, b)
                cormat <- cov2rmaslope(lowerdiag2mat(sigma_samples[i,ab], col_names = param_names,
                                                     hier = TRUE))
                cormat_flat <- flatten_matrix(cormat, diag = FALSE)
                bb <- j * nm
                aa <- bb - nm + 1
                aabb <- seq(aa, bb)
                cormat_samples[i,aabb] <- cormat_flat
            }
        } else {
            # Simple multivariate model -- nparam
            cormat <- cov2rmaslope(lowerdiag2mat(sigma_samples[i,], col_names = param_names))
            cormat_flat <- flatten_matrix(cormat, diag = FALSE)
            cormat_samples[i,] <- cormat_flat
        }
    }
    dn <- dimnames(sigma_samples)
    dn[[2]] <- corr_names
    dimnames(cormat_samples) <- dn
    return(cormat_samples)
}

#' Calculate reduced major axis regression slopes from a covariance matrix
#'
#' @param m Variance-covariance matrix
#' @return
#' @author Alexey Shiklomanov
#' @examples
#' cov2rmaslope(cov(iris[, -5]))
#' @export
cov2rmaslope <- function(m) {
    # RMA slope = sign(cov(xy)) * sd(y) / sd(x)
    s <- sqrt(diag(m))
    ones <- matrix(1, nrow(m), ncol(m))
    diag(ones) <- s
    a <- ones %*% diag(s)
    b <- diag(s) %*% ones
    out <- sign(m) * b / a
    out
}
