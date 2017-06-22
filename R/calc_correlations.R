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
    cormat_samples <- sigma_samples
    if (length(dims) == 4) {
        for (i in seq_len(dims[1])) {
            if (length(dims) == 4) {
            # Hierarchical model -- nsamp, ngroup, nparam, nparam
                for (j in seq_len(dims[2])) {
                    cormat_samples[i,j,,] <- cov2cor(sigma_samples[i,j,,])
                }
            } else {
            # Simple multivariate model -- nsamp, nparam, nparam
                cormat_samples[i,,] <- cov2cor(sigma_samples[i,,])
            }
        }
    }
    return(cormat_samples)
}
