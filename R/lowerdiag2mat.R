#' Convert a vector of the lower diagonoal of a symmetrical matrix to a matrix
#'
#' Storing just the lower diagonal is an efficient way to store MCMC samples of
#' a matrix (and, in fact, is how matrices are stored by the samplers in
#' [fit_mvnorm()] and [fit_mvnorm_hier()]).
#'
#' @param vec Vector of lower diagonal values
#' @param col_names (Logical) If `TRUE` (default), extract row and column names from
#'   vector names formatted as `"row..column"`.
#' @param corr (Logical) If `TRUE`, assume this is a correlation matrix where
#'   the diagonal is fixed at 1 and therefore not stored.
#' @param colorder Optional numeric or character vector specifying the desired
#'   column order.
#' @param hier (Logical) Whether the vector names also include a group name.
#'   Only used if `col_names` is `TRUE`.
#' @export
lowerdiag2mat <- function(vec, col_names = TRUE, corr = FALSE, colorder = NULL, hier = FALSE) {
  # Determine matrix dimensions
  nvec <- length(vec)
  nmat <- 0.5 * (sqrt(8 * nvec + 1) - 1)
  stopifnot(nmat %% 1 == 0)
  if (corr) nmat <- nmat + 1

  # Populate matrix
  mat <- diag(nmat)
  mat[lower.tri(mat, diag = !corr)] <- vec
  mat <- mat + t(mat) - diag(diag(mat))

  # Try to figure out matrix dimnames from vector names
  if (isTRUE(col_names)) {
    splitnames <- do.call(rbind, strsplit(names(vec), split = varsep_esc))
    if (hier) {
      group <- unique(splitnames[,1])
      stopifnot(length(group) == 1)
      pars <- splitnames[seq_len(nmat), 2]
      col_names <- paste(group, pars, sep = varsep)
    } else {
      col_names <- splitnames[seq_len(nmat), 1]
    }
    stopifnot(length(col_names) == nmat)
    rownames(mat) <- colnames(mat) <- col_names
  }
  if (!is.null(colorder)){
    mat <- mat[colorder, colorder]
  }
  return(mat)
}
