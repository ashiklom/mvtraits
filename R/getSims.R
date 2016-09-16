# Functions to retrieve outputs from RData files

#' @export
get_sims <- function(filename, vars){
    library(rstan)
    load(filename)
    stopifnot(length(vars) == 1)
    out_sims <- rstan::extract(out, vars)[[1]]
    dimnames(out_sims)[[2]] <- traits_nolog
    remove(out)
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

#' @export
selectFromMatrix <- function(mat, variable) {
    mat[, grep(variable, colnames(mat))]
}

#' @export
getPattern <- function(pattern, string) gsub(pattern, "\\1", string)

#' @export
assignTraitNames <- function(mat){
    oldnames <- colnames(mat)
    newnames <- sapply(oldnames, replaceName)
    return(newnames)
}

#' @export
replaceName <- function(colname){
    prefix <- getPattern("(.*)\\[.*", colname)
    if (grepl("mu_trait|sigma2_obvs", colname)) {
        trait_index <- as.numeric(getPattern(".*\\[([[:digit:]]+)\\]", colname))
        newname <- paste(prefix, traits[trait_index], sep=".")
    } else if (grepl("Sigma_trait", colname)) {
        trait1 <- as.numeric(getPattern(".*\\[([[:digit:]]+),[[:digit:]]+\\]", colname))
        trait2 <- as.numeric(getPattern(".*\\[[[:digit:]]+,([[:digit:]]+)\\]", colname))
        newname <- paste(prefix, traits[trait1], traits[trait2], sep=".")
    } else if (grepl("mu_pft_trait", colname)) {
        pft_number <- as.numeric(getPattern(".*\\[([[:digit:]]+),[[:digit:]]+\\]", colname))
        trait_index <- as.numeric(getPattern(".*\\[[[:digit:]]+,([[:digit:]]+)\\]", colname))
        newname <- paste(prefix, pft.names[pft_number], traits[trait_index], sep=".")
    } else if (grepl("Sigma_pft", colname)) {
        pft_number <- as.numeric(getPattern(".*\\[([[:digit:]]+),[[:digit:]]+,[[:digit:]]+\\]", colname))
        trait1 <- as.numeric(getPattern(".*\\[[[:digit:]]+,([[:digit:]]+),[[:digit:]]+\\]", colname))
        trait2 <- as.numeric(getPattern(".*\\[[[:digit:]]+,[[:digit:]]+,([[:digit:]]+)\\]", colname))
        newname <- paste(prefix, pft.names[pft_number], traits[trait1], traits[trait2], sep=".")
    }
    return(newname)
}


