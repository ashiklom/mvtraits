#' @export
gibbs_default_priors <- function(nparam, ngroup = NULL) {
    prior <- list(mu_global = rep(0, nparam), 
                  Sigma_global = diag(1.0, nparam),
                  v_global = 0,
                  S_global = diag(1.0, nparam))
    if (!is.null(ngroup) && !is.na(ngroup)) {
        group <- list(mu_group = matrix(0, ngroup, nparam), 
                      Sigma_group = matrep(diag(1.0, nparam), ngroup),
                      v_group = rep(0, ngroup),
                      S_group = matrep(diag(1.0, nparam)))
        prior <- c(prior, group)
    }
    return(prior)
}

#' @export
Wishart_prior_param <- function(mean, sd, nparam) {
    v0 <- max(0, 2 * mean ^ 2 / sd ^ 2 - nparam - 1)
    df <- v0 + nparam + 1
    s <- mean / df
    S0 <- diag(s, nparam)
    return(list(v0 = v0, S0 = S0))
}

#' @export
matrep <- function(mat, n) {
    aperm(replicate(ngroup, mat), c(3, 1, 2))
}
