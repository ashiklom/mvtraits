#' @export
chain2matrix_hier <- function(chain) {
    mu_global_mat <- chain[['mu_global']]
    colnames(mu_global_mat) <- paste('mu', 'global', colnames(mu_global_mat), sep = varsep)
    sigma_global_mat <- chain[['Sigma_global']]
    colnames(sigma_global_mat) <- paste('Sigma', 'global', colnames(sigma_global_mat), sep = varsep)
    mu_group_mat <- chain[['mu_group']]
    colnames(mu_group_mat) <- paste('mu', colnames(mu_group_mat), sep = varsep)
    sigma_group_mat <- chain[['Sigma_group']]
    colnames(sigma_group_mat) <- paste('Sigma', colnames(sigma_group_mat), sep = varsep)
    if ('Corr_global' %in% names(chain) && 'Corr_group' %in% names(chain)) {
        corr_global_mat <- chain[['Corr_global']]
        colnames(corr_global_mat) <- paste('Corr', 'global', colnames(corr_global_mat), sep = varsep)
        corr_group_mat <- chain[['Corr_group']]
        colnames(corr_group_mat) <- paste('Corr', colnames(corr_group_mat), sep = varsep)
        params_mat <- cbind(mu_global_mat, sigma_global_mat, corr_global_mat, 
                            mu_group_mat, sigma_group_mat, corr_group_mat)
    } else {
        params_mat <- cbind(mu_global_mat, sigma_global_mat, mu_group_mat, sigma_group_mat)
    }
    return(params_mat)
}
