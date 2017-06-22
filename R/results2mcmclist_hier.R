#' @export
chain2matrix_hier <- function(chain) {
    mu_global_mat <- chain[['mu_global']]
    colnames(mu_global_mat) <- paste('mu', 'global', colnames(mu_global_mat), sep = varsep)
    sigma_global_mat <- flatten_sigma(chain[['Sigma_global']])
    colnames(sigma_global_mat) <- paste('Sigma', 'global', colnames(sigma_global_mat), sep = varsep)
    mu_group_mat <- flatten_mu_group(chain[['mu_group']])
    colnames(mu_group_mat) <- paste('mu', colnames(mu_group_mat), sep = varsep)
    sigma_group_mat <- flatten_sigma_group(chain[['Sigma_group']])
    colnames(sigma_group_mat) <- paste('Sigma', colnames(sigma_group_mat), sep = varsep)
    params_mat <- cbind(mu_global_mat, sigma_global_mat, mu_group_mat, sigma_group_mat)
    return(params_mat)
}

flatten_mu_group <- function(mu_group_samples) {
    dims <- dim(mu_group_samples)
    dimsnames <- dimnames(mu_group_samples)
    nsamp <- dims[1]
    ngroup <- dims[2]
    nparam <- dims[3]
    ncolumn <- ngroup * nparam
    group_names <- dimsnames[[2]]
    param_names <- dimsnames[[3]]
    mu_group_mat_names <- paste(rep(group_names, nparam), rep(param_names, each = ngroup), sep = varsep)
    mu_group_mat <- matrix(mu_group_samples, nrow = nsamp, ncol = ncolumn)
    colnames(mu_group_mat) <- mu_group_mat_names
    return(mu_group_mat)
}

flatten_sigma_group <- function(sigma_group_samples) {
    ngroup <- dim(sigma_group_samples)[2]
    dimsnames <- dimnames(sigma_group_samples)
    group_names <- dimsnames[[2]]
    param_names <- dimsnames[[3]]
    group_list <- lapply(seq_len(ngroup), function(i) flatten_sigma(sigma_group_samples[,i,,]))
    ncolumn <- ncol(group_list[[1]])
    group_mat <- do.call(cbind, group_list)
    colnames(group_mat) <- paste(rep(group_names, each = ncolumn), colnames(group_mat), sep = varsep)
    return(group_mat)
}
