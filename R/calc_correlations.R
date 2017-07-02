#' @export
add_correlations <- function(results_list, hier = FALSE, ngroups = NULL) {
    vars <- names(results_list[[1]])
    for (i in seq_along(results_list)) {
        if (hier) {
            results_list[[i]][['Corr_global']] <- calc_correlations(results_list[[i]][['Sigma_global']])
            results_list[[i]][['Corr_group']] <- calc_correlations(results_list[[i]][['Sigma_group']], hier = TRUE, ngroups = ngroups)
        } else {
            results_list[[i]][['Corr']] <- calc_correlations(results_list[[i]][['Sigma']])
        }
    }
    return(results_list)
}

#' @export
calc_correlations <- function(sigma_samples, hier = FALSE, ngroups = NULL) {
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
    for (i in seq_len(dims[1])) {
        if (hier) {
            # Hierarchical model -- ngroup (slow) x nparam (fast)
            cormat_names <- character(cdims[2])
            for (j in seq_len(ngroups)) {
                b <- j * nvec
                a <- b - nvec + 1
                ab <- seq(a, b)
                cormat <- cov2cor(lowerdiag2mat(sigma_samples[i,ab], hier = TRUE))
                cormat_flat <- flatten_matrix(cormat, diag = FALSE)
                bb <- j * nm
                aa <- bb - nm + 1
                aabb <- seq(aa, bb)
                cormat_samples[i,aabb] <- cormat_flat
                cormat_names[aabb] <- names(cormat_flat)
            }
        } else {
            # Simple multivariate model -- nparam
            cormat <- cov2cor(lowerdiag2mat(sigma_samples[i,]))
            cormat_flat <- flatten_matrix(cormat, diag = FALSE)
            cormat_samples[i,] <- cormat_flat
        }
    }
    dn <- dimnames(sigma_samples)
    dn[[2]] <- if (hier) cormat_names else names(cormat_flat)
    dimnames(cormat_samples) <- dn
    return(cormat_samples)
}
