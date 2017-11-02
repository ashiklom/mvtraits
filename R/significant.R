#' Test if pairwise correlations are significant
#'
#' @param matlist List of covariance (or correlation) matrices, such as that
#' returned by [fit_mvnorm()] or [fit_mvnorm_hier()] summary functions.
#' @param xvar Name of x variable as string (must match row name in `matlist`)
#' @param yvar Name of y variable as string (must match column name in
#' `matlist`)
#' @param lo Name of lower significance statistic cutoff (default = "2.5\%")
#' @param hi Name of upper significance statistic cutoff (default = "97.5\%")
test_significant <- function(matlist, xvar, yvar, lo = "2.5%", hi = "97.5%") {
  xy <- c(xvar, yvar)
  m1 <- matlist[[lo]][xy, xy]
  m2 <- matlist[[hi]][xy, xy]
  all(sign(m1) == sign(m2))
}
