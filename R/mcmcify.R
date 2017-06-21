#' @export
results2mcmclist_multi <- function(results_list) {
    matrified <- lapply(results_list, chain2matrix_multi)
    matrix_mcmc <- lapply(matrified, coda::as.mcmc)
    mcmc_list <- coda::as.mcmc.list(matrix_mcmc)
    return(mcmc_list)
} 

chain2matrix_multi <- function(chain) {
    mu_mat <- chain[['mu']]
    colnames(mu_mat) <- paste('mu', colnames(mu_mat), sep = '.')
    sigma_mat <- flatten_sigma(chain[['Sigma']])
    colnames(sigma_mat) <- paste('Sigma', colnames(sigma_mat), sep = '.')
    params_mat <- cbind(mu_mat, sigma_mat)
    return(params_mat)
}

flatten_sigma <- function(sigma_samples) {
    dims <- dim(sigma_samples)
    nparam <- tail(dims, 1)
    param_names <- dimnames(sigma_samples)[[length(dims)]]
    tri_inds <- which(upper.tri(diag(nparam), diag = TRUE), arr.ind = TRUE)
    sigma_mat_names <- paste(param_names[tri_inds[,'row']], param_names[tri_inds[,'col']], sep = '.')
    sigma_mat <- t(apply(sigma_samples, 1, function(x) x[upper.tri(x, diag = TRUE)]))
    colnames(sigma_mat) <- sigma_mat_names
    return(sigma_mat)
}

chain2matrix_hier <- function(chain) {
    mu_global_mat <- chain[['mu']]
    colnames(mu_global_mat) <- paste('mu_global', colnames(mu_global_mat), sep = '.')
    sigma_global_mat <- flatten_sigma(chain[['Sigma_global']])
    colnames(sigma_global_mat) <- paste('Sigma_global', colnames(sigma_global_mat), sep = '.')
    params_mat <- cbind(mu_mat, sigma_mat)
    return(params_mat)
}

