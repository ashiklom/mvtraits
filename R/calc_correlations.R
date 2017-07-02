#' @export
add_correlations <- function(results_list) {
    vars <- names(results_list[[1]])
    hier <- 'Sigma_global' %in% vars && 'Sigma_group' %in% vars
    for (i in seq_along(results_list)) {
        if (hier) {
            results_list[[i]][['Corr_global']] <- calc_correlations(results_list[[i]][['Sigma_global']])
            results_list[[i]][['Corr_group']] <- calc_correlations(results_list[[i]][['Sigma_group']])
        } else {
            results_list[[i]][['Corr']] <- calc_correlations(results_list[[i]][['Sigma']])
        }
    }
    return(results_list)
}

#' @export
calc_correlations <- function(sigma_samples) {
    dims <- dim(sigma_samples)
    if (length(dims) == 4) {
        hier <- TRUE
    } else {
        hier <- FALSE
    }
    cdims <- dims
    nvec <- dims[2]
    nc <- 0.5 * (sqrt(8 * nvec + 1) - 1)
    nc2 <- nc - 1
    cdims[2] <- nc2 * (nc2 + 1) / 2
    cormat_samples <- array(0, cdims)
    for (i in seq_len(dims[1])) {
        if (hier) {
            # Hierarchical model -- nsamp, ngroup, nparam, nparam
            for (j in seq_len(dims[2])) {
                cormat_samples[i,j,,] <- cov2cor(sigma_samples[i,j,,])
            }
        } else {
            # Simple multivariate model -- nsamp, nparam, nparam
            cormat <- cov2cor(lowerdiag2mat(sigma_samples[i,]))
            cormat_flat <- flatten_matrix(cormat, diag = FALSE)
            cormat_samples[i,] <- cormat_flat
        }
    }
    dn <- dimnames(sigma_samples)
    dn[[2]] <- names(cormat_flat)
    dimnames(cormat_samples) <- dn
    return(cormat_samples)
}
