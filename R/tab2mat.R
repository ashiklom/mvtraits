#' @export
tab2mat <- function(dat, colname = "Mean") {
    library(data.table)
    stopifnot(is.data.table(dat))
    vec <- dat[, colname, with = FALSE]
    mat <- lowerdiag2mat(vec)
    return(mat)
}

#' @export
lowerdiag2mat <- function(vec) {
    nvec <- length(vec)
    stopifnot(nvec %% 1 == 0)
    nmat <- 0.5 * (sqrt(8 * nvec + 1) - 1)
    stopifnot(nmat %% 1 == 0)
    mat <- matrix(numeric(), nmat, nmat)
    mat[lower.tri(mat, diag=TRUE)] <- vec
    mat[upper.tri(mat)] <- mat[lower.tri(mat)]
    return(mat)
}
