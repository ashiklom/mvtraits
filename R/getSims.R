# Functions to retrieve outputs from RData files
selectFromMatrix <- function(mat, variable) mat[, grep(variable, colnames(mat))]

getPattern <- function(pattern, string) gsub(pattern, "\\1", string)

replaceName <- function(colname){
    if (grepl("mu_.*trait|sigma2_obvs", colname)) {
        prefix <- getPattern("(.*)\\[.*", colname)
        trait_index <- as.numeric(getPattern(".*\\[([[:digit:]]+)\\]", colname))
        newname <- paste(prefix, traits[trait_index], sep=".")
    } else if (grepl("Sigma", colname)) {
        prefix <- getPattern("(.*)\\[.*", colname)
        trait1 <- as.numeric(getPattern(".*\\[([[:digit:]]+),[[:digit:]]+\\]", colname))
        trait2 <- as.numeric(getPattern(".*\\[[[:digit:]]+,([[:digit:]]+)\\]", colname))
        newname <- paste(prefix, traits[trait1], traits[trait2], sep=".")
    }
    return(newname)
}

assignTraitNames <- function(mat){
    oldnames <- colnames(mat)
    newnames <- sapply(oldnames, replaceName)
    return(newnames)
}

get_sims <- function(filename, vars){
    library(runjags)
    library(dplyr)
    load(filename)
    matchvars <- paste(vars, collapse="|")
    out_sims <- out$mcmc %>% as.matrix() %>% selectFromMatrix(vars)
    colnames(out_sims) <- assignTraitNames(out_sims)
    remove(out)
    return(out_sims)
}

get_sims_list <- function(filenames, vars){
    l <- list()
    for(f in filenames){
        pft_number <- gsub(".*_([[:digit:]]+).*", "\\1", f)
        pft_name <- pft.names[pft_number]
        l[[pft_name]] <- get_sims(f, vars)
    }
    return(l)
}
