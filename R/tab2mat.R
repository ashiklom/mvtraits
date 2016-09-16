#' @export
tab2mat <- function(dat, colname = "Mean", ...) {
    library(data.table)
    stopifnot(is.data.table(dat))
    vec <- dat[[colname]]
    mat <- lowerdiag2mat(vec, ...)
    return(mat)
}

#' @export
lowerdiag2mat <- function(vec, corr = FALSE) {
    nvec <- length(vec)
    stopifnot(nvec %% 1 == 0)
    nmat <- 0.5 * (sqrt(8 * nvec + 1) - 1)
    stopifnot(nmat %% 1 == 0)
    if (corr) nmat <- nmat + 1
    mat <- matrix(numeric(), nmat, nmat)
    mat[lower.tri(mat, diag=!corr)] <- vec
    mat[upper.tri(mat)] <- mat[lower.tri(mat)]
    if (corr) diag(mat) <- 1
    return(mat)
}
