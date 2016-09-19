# Functions to retrieve outputs from RData files

#' @export
get_sims <- function(filename, vars, samples = 5000){
    library(rstan)
    load(filename)
    ntrait <- length(traits_nolog)
    nt <- 1:ntrait
    out_sims <- rstan::extract(out, vars)
    s <- function(x) sample.int(dim(x)[1], samples)
    remove(out)
    if ("mu" %in% vars) {
        n <- s(out_sims[['mu']])
        out_sims[["mu"]] <- out_sims[["mu"]][n,nt]
        colnames(out_sims[["mu"]]) <- traits_nolog
    }
    if ("mu_global" %in% vars) {
        n <- s(out_sims[['mu_global']])
        out_sims[["mu_global"]] <- out_sims[["mu_global"]][n,nt]
        colnames(out_sims[["mu_global"]]) <- traits_nolog
    }
    if ("mu_pft" %in% vars) {
        n <- s(out_sims[['mu_pft']])
        out_sims[["mu_pft"]] <- out_sims[["mu_pft"]][n,,nt]
        dimnames(out_sims[["mu_pft"]])[[2]] <- pft.names
        dimnames(out_sims[["mu_pft"]])[[3]] <- traits_nolog
    }
    if ("sigma" %in% vars) {
        n <- s(out_sims[['sigma2']])
        out_sims[["sigma2"]] <- out_sims[["sigma2"]][n,nt,nt]
        dimnames(out_sims[["sigma2"]])[[2]] <- traits_nolog
    }
    if ("Sigma" %in% vars) {
        n <- s(out_sims[['Sigma']])
        out_sims[["Sigma"]] <- out_sims[["Sigma"]][n,nt,nt]
        dimnames(out_sims[["Sigma"]])[[2]] <- traits_nolog
        dimnames(out_sims[["Sigma"]])[[3]] <- traits_nolog
    }
    if ("Sigma_global" %in% vars) {
        n <- s(out_sims[['Sigma_global']])
        out_sims[["Sigma_global"]] <- out_sims[["Sigma_global"]][n,nt,nt]
        dimnames(out_sims[["Sigma_global"]])[[2]] <- traits_nolog
        dimnames(out_sims[["Sigma_global"]])[[3]] <- traits_nolog
    }
    if ("Sigma_pft" %in% vars) {
        n <- s(out_sims[['Sigma_pft']])
        out_sims[["Sigma_pft"]] <- out_sims[["Sigma_pft"]][n,,nt,nt]
        dimnames(out_sims[["Sigma_pft"]])[[2]] <- pft.names
        dimnames(out_sims[["Sigma_pft"]])[[3]] <- traits_nolog
        dimnames(out_sims[["Sigma_pft"]])[[4]] <- traits_nolog
    }
    return(out_sims)
}

#' @export
get_sims_list <- function(filenames, vars){
    l <- list()
    for(f in filenames){
        pft_number <- gsub(".*_([[:digit:]]+).*", "\\1", f) %>%
            as.numeric
        pft_name <- pft.names[pft_number]
        l[[pft_name]] <- get_sims(f, vars)
    }
    return(l)
}
