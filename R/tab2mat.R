#' @export
tab2mat <- function(dat, colname = "Mean", ...) {
    library(data.table)
    stopifnot(is.data.table(dat))
    vec <- dat[[colname]]
    names(vec) <- dat[["trait"]]
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
    namesvec <- names(vec)
    splitnames <- strsplit(namesvec, split="\\.")
    allpars <- unique(unlist(splitnames))
    stopifnot(length(allpars) == nmat)
    mat <- matrix(numeric(), nmat, nmat)
    rownames(mat) <- colnames(mat) <- allpars
    for (i in seq_along(splitnames)) {
        mat[splitnames[[i]][1],splitnames[[i]][2]] <- vec[namesvec[i]]
        mat[splitnames[[i]][2],splitnames[[i]][1]] <- vec[namesvec[i]]
    }
    if (corr) diag(mat) <- 1
    return(mat)
}
