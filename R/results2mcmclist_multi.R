#' @export
results2mcmclist <- function(results_list, chainfun) {
    matrified <- lapply(results_list, chainfun)
    matrix_mcmc <- lapply(matrified, coda::as.mcmc)
    mcmc_list <- coda::as.mcmc.list(matrix_mcmc)
    return(mcmc_list)
} 

#' @export
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
    tri_inds <- which(lower.tri(diag(nparam), diag = TRUE), arr.ind = TRUE)
    sigma_mat_names <- paste(param_names[tri_inds[,'row']], param_names[tri_inds[,'col']], sep = '.')
    sigma_mat <- t(apply(sigma_samples, 1, function(x) x[lower.tri(x, diag = TRUE)]))
    colnames(sigma_mat) <- sigma_mat_names
    return(sigma_mat)
}

