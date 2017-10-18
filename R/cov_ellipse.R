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
  # calculate axes
  radius <- sqrt(max_val)
  x0 <- mean[1]
  y0 <- mean[2]
  x1 <- x0 - radius * cos(angle)
  x2 <- x0 + radius * cos(angle)
  y1 <- y0 - radius * sin(angle)
  y2 <- y0 + radius * sin(angle)

  tibble::tibble(
    angle = angle,
    radius = sqrt(max_val),
    center_x = x0,
    center_y = y0,
    x1 = x1,
    x2 = x2,
    y1 = y1,
    y2 = y2,
    ellipse_x = list(r_ellipse[, 1] + mean[1]),
    ellipse_y = list(r_ellipse[, 2] + mean[2])
  )
}

#' Draw ellipse
#'
#' @param ellipse_dat Ellipse data.frame, as returned by [ellipse_axes]
#' @inheritParams graphics::lines
#' @param ... Additional arguments to [graphics::plot]
#' @export
draw_ellipse <- function(ellipse_dat, col = 1, lwd = 1, lty = 1, ...) {
  x <- ellipse_dat[[1, "ellipse_x"]]
  y <- ellipse_dat[[1, "ellipse_y"]]
  if (dev.cur() == 1) {
    # Plot has not been called. Create a new one.
    plot(x, y, type = "n", ...)
  }
  lines(x, y, col = col, lwd = lwd, lty = lty)
}

#' Draw ellipse major axis
#'
#' @inheritParams draw_ellipse
#' @inheritParams graphics::lines
#' @inheritParams graphics::plot
#' @export
draw_majoraxis <- function(ellipse_dat, col = 1, lwd = 1, lty = 1, pch = 19, ...) {
  with(ellipse_dat, {
    if (dev.cur() == 1) {
      plot(c(x1, x2), c(y1, y2), type = "n", ...)
    }
    segments(x1, y1, x2, y2, col = col, lwd = lwd, lty = lty)
    points(center_x, center_y, pch = pch)
  })
}
