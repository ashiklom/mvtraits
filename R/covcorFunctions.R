vecToCovMat <- function(vec) {
    cols <- sqrt(length(vec))
    if (cols %% 1 != 0) stop("Not a square matrix")
    covmat <- matrix(vec, cols)
    return(covmat)
}

covToCor <- function(vec){
    covmat <- vecToCovMat(vec)
    cormat <- cov2cor(covmat)
    corvec <- c(cormat)
    return(corvec)
}

covToCor.global <- function(mat){
    cormat <- t(apply(mat, 1, covToCor))
    colnames(cormat) <- colnames(mat)
    return(cormat)
}

