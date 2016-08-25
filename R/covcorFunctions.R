covToCor <- function(vec){
    cols <- sqrt(length(vec))
    if (cols %% 1 != 0) stop("Not a square matrix")
    covmat <- matrix(vec, cols)
    cormat <- cov2cor(covmat)
    corvec <- c(covmat)
    return(corvec)
}

covToCor.global <- function(mat){
    cormat <- t(apply(cov.global, 1, covToCor))
    colnames(cormat) <- colnames(mat)
    return(cormat)
}
