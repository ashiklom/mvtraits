#' @export
results2mcmclist <- function(results_list, type) {
    if (type == 'multi') {
        chainfun <- chain2matrix_multi
    } else if (type == 'hier') {
        chainfun <- chain2matrix_hier
    } else {
        stop('Unknown type ', type)
    }
    matrified <- lapply(results_list, chainfun)
    matrix_mcmc <- lapply(matrified, coda::as.mcmc)
    mcmc_list <- coda::as.mcmc.list(matrix_mcmc)
    return(mcmc_list)
} 

#' @export
chain2matrix_multi <- function(chain) {
    mu_mat <- chain[['mu']]
    colnames(mu_mat) <- paste('mu', colnames(mu_mat), sep = varsep)
    sigma_mat <- chain[['Sigma']]
    colnames(sigma_mat) <- paste('Sigma', colnames(sigma_mat), sep = varsep)
    if ('Corr' %in% names(chain)) {
        corr_mat <- chain[['Corr']]
        colnames(corr_mat) <- paste('Corr', colnames(corr_mat), sep = varsep)
        params_mat <- cbind(mu_mat, sigma_mat, corr_mat)
    } else {
        params_mat <- cbind(mu_mat, sigma_mat)
    }
    return(params_mat)
}

lowertri_names <- function(col_names, diag = TRUE) {
    n <- length(col_names)
    tri_inds <- which(lower.tri(diag(n), diag = diag), arr.ind = TRUE)
    vec_names <- paste(col_names[tri_inds[,'row']], col_names[tri_inds[,'col']], sep = varsep)
    return(vec_names)
}

flatten_matrix <- function(mat, diag = TRUE) {
    stopifnot(nrow(mat) == ncol(mat))
    vec <- mat[lower.tri(mat, diag = diag)]
    vnames <- colnames(mat)
    if (!is.null(vnames)) {
        names(vec) <- lowertri_names(vnames, diag = diag)
    }
    return(vec)
}

flatten_sigma <- function(sigma_samples, diag = TRUE) {
    dims <- dim(sigma_samples)
    nparam <- tail(dims, 1)
    param_names <- dimnames(sigma_samples)[[length(dims)]]
    tri_inds <- which(lower.tri(diag(nparam), diag = diag), arr.ind = TRUE)
    sigma_mat_names <- paste(param_names[tri_inds[,'row']], param_names[tri_inds[,'col']], sep = varsep)
    sigma_mat <- t(apply(sigma_samples, 1, function(x) x[lower.tri(x, diag = diag)]))
    colnames(sigma_mat) <- sigma_mat_names
    return(sigma_mat)
}

varsep <- '..'
varsep_esc <- gsub('(\\W)', '\\\\\\1', varsep, perl = TRUE)
