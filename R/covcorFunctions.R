#' @export
vecToMat <- function(vec) {
    cols <- sqrt(length(vec))
    if (cols %% 1 != 0) stop("Not a square matrix")
    covmat <- matrix(vec, cols)
    return(covmat)
}

#' @export
covToCor <- function(vec){
    covmat <- vecToMat(vec)
    cormat <- cov2cor(covmat)
    corvec <- c(cormat)
    return(corvec)
}

#' @export
covToCor.global <- function(mat){
    cormat <- t(apply(mat, 1, covToCor))
    colnames(cormat) <- colnames(mat)
    return(cormat)
}
