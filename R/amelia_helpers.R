#' @export
arraycorr <- function(arr) {
    out <- arr
    for (i in seq_len(dim(arr)[3])) {
        out[,,i] <- cov2cor(arr[,,i])
    }
    return(out)
}

#' @export
vec2df <- function(vec, ...) tibble(index = names(vec), value = vec, ...)

#' @export
mat2df <- function(mat, diag = TRUE, ...) {
    pars <- rownames(mat)
    lt <- which(lower.tri(mat, diag = diag), arr.ind = TRUE)
    matvec <- mat[lt]
    names(matvec) <- paste(pars[lt[,1]], pars[lt[,2]], sep = '..')
    vec2df(matvec, ...)
}

#' @export
imp2 <- function(dat, trait_order) {
    dat %>% 
        filter(sapply(imputed, class) == 'amelia') %>% 
        mutate(mu_all = map(imputed, 'mu'),
               mu_all = map2(mu_all, data, ~"rownames<-"(.x, colnames(.y))),
               mu_all = map(mu_all, ~.x[trait_order,]),
               mu_means = map(mu_all, rowMeans),
               Sigma_all = map(imputed, 'covMatrices'),
               Sigma_all = map2(Sigma_all, data, ~"dimnames<-"(.x, list(colnames(.y), colnames(.y), NULL))),
               Sigma_all = map(Sigma_all, ~.x[trait_order,trait_order,]),
               Sigma_means = map(Sigma_all, apply, 1:2, mean),
               Corr_all = map(Sigma_all, arraycorr),
               Corr_means = map(Corr_all, apply, 1:2, mean))
}

#' @export
proc <- function(imp2) {
    imp2 %>% 
        mutate(mu_df = map(mu_means, vec2df, variable = 'mu'),
               sigma_df = map(Sigma_means, mat2df, variable = 'sigma'),
               corr_df = map(Corr_means, mat2df, variable = 'corr', diag = FALSE)) %>% 
        mutate(alldat = pmap(list(mu_df, sigma_df, corr_df), bind_rows)) %>% 
        select(pft, alldat) %>% 
        unnest()
}

#' @export
means_wide <- function(proc, trait_order) {
    means_wide_all <- proc %>% 
        filter(variable == 'mu') %>% 
        select(-variable) %>% 
        mutate(index = factor(index, trait_order)) %>% 
        spread(index, value)
    means_wide <- select(means_wide_all, -pft)
    return(means_wide)
}

#' @export
bind_all <- function(means_wide, proc, trait_order) {
    glob_mean <- vec2df(colMeans(means_wide), pft = 'global', variable = 'mu')
    glob_cov <- mat2df(cov(means_wide), pft = 'global', variable = 'sigma')
    glob_cor <- mat2df(cor(means_wide), diag = FALSE, pft = 'global', variable = 'corr')
    glob <- bind_rows(glob_mean, glob_cov, glob_cor)
    glob_all <- full_join(glob, proc) %>% 
        mutate(pft = factor(pft, levels = c('global', levels(proc[['pft']])))) %>% 
        separate(index, c('xvar', 'yvar'), sep = '\\.\\.', remove = FALSE) %>% 
        mutate(xvar = factor(xvar, levels = trait_order),
               yvar = factor(yvar, levels = trait_order))
    return(glob_all)
}
