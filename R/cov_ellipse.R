#' Calculate axes for covariance matrix ellipse
#'
#' @param cov Variance-covariance matrix
#' @param mean Mean vector
#' @param prob Ellipse probability quantile
#' @export
ellipse_axes <- function(cov, mean = c(0, 0), prob = 0.95) {
  stopifnot(
    length(mean) == ncol(cov),
    length(mean) == 2,
    is.matrix(cov)
  )
  eig <- eigen(cov)
  max_val <- max(eig$values)
  imax <- which(eig$values == max_val)
  max_vec <- eig$vectors[, imax]
  min_val <- min(eig$values)
  imin <- which(eig$values == min_val)
  min_vec <- eig$vectors[, imin]
  angle <- atan2(max_vec[2], max_vec[1])
  theta <- seq(0, 2 * pi, 0.01)
  sqrt_chisq <- sqrt(qchisq(prob, 2))
  a <- sqrt_chisq * sqrt(max_val)
  b <- sqrt_chisq * sqrt(min_val)
  ellipse_x <- a * cos(theta)
  ellipse_y <- b * sin(theta)
  # rotation matrix
  rmat <- cbind(c(cos(angle), -sin(angle)), c(sin(angle), cos(angle)))
  # apply rotation, and add means to each row
  r_ellipse <- cbind(ellipse_x, ellipse_y) %*% rmat
  tibble::tibble(
    angle = angle,
    radius = sqrt(max_val),
    ellipse_x = list(r_ellipse[, 1] + mean[1]),
    ellipse_y = list(r_ellipse[, 2] + mean[2])
  )
}
