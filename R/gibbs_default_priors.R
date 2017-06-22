#' @export
gibbs_default_priors <- function(nparam, ngroup = NULL) {
    prior <- list(mu_global = rep(0, nparam), 
                  Sigma_global = diag(1.0, nparam),
                  v_global = 0,
                  S_global = diag(1.0, nparam))
    if (!is.null(ngroup) && !is.na(ngroup)) {
        group <- list(mu_group = matrix(0, ngroup, nparam), 
                      Sigma_group = aperm(replicate(ngroup, diag(1.0, nparam)), c(3, 1, 2)),
                      v_group = rep(0, ngroup),
                      S_group = aperm(replicate(ngroup, diag(1.0, nparam)), c(3, 1, 2)))
        prior <- c(prior, group)
    }
    return(prior)
}
