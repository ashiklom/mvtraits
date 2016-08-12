# Functions to retrieve outputs from RData files
selectFromMatrix <- function(mat, variable) mat[, grep(variable, colnames(mat))]

get_sims <- function(filename, vars){
    library(runjags)
    library(dplyr)
    load(filename)
    matchvars <- paste(vars, collapse="|")
    out_sims <- out$mcmc %>% as.matrix() %>% selectFromMatrix(vars)
    remove(out)
    return(out_sims)
}

get_sims_list <- function(filenames, vars){
    l <- list()
    for(f in seq_along(filenames)){
        l[[f]] <- get_sims(filenames[f], vars)
    }
    names(l) <- gsub("(.*)\\.Rdata", "\\1", filenames) 
    return(l)
}
