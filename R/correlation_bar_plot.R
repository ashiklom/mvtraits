#' Draw matrix of correlation bar plots
#'
#' @param dat_lower Data for lower triangle
#' @param dat_upper Data for upper triangle
#' @param param_order Character vector of parameter values indicating order
#' @param ... Additional arguments to [single_bar_plot]
bar_plot_matrix <- function(dat_lower, dat_upper,
                            param_order = NULL, significant_only = TRUE, ...) {
  if (is.null(param_order)) {
    param_order <- union(dat_lower[["xvar"]], dat_lower[["yvar"]])
  }
  nparam <- length(param_order)
  nr <- nrow(dat_lower)
  nrseq <- seq_len(nr)
  xind_lower <- xind_upper <- match(dat_lower[["xvar"]], param_order)
  yind_lower <- yind_upper <- match(dat_lower[["yvar"]], param_order)
  dat_lower_r <- dat_lower[order(yind_lower, xind_lower), ]
  dat_upper_r <- dat_upper[order(yind_upper, xind_upper), ]

  lay_mat <- matrix(0, nparam, nparam)
  lay_mat[lower.tri(lay_mat)] <- nrseq
  lay_mat[upper.tri(lay_mat)] <- nr + nrseq
  diag(lay_mat) <- 2 * nr + nrseq
  layout(lay_mat)
  lapply(
    dat_lower_r[["data"]],
    single_bar_plot,
    significant_only = significant_only,
    ...
  )
  lapply(
    dat_upper_r[["data"]],
    single_bar_plot,
    significant_only = significant_only,
    ...
  )
  lapply(param_order, text_plot, cex = 2)
}

#' Single bar plot panel
#' 
#' @param dat Data for panel
#' @param significant_only If `TRUE`, hide columns that are not significantly 
#' difference from zero.
single_bar_plot <- function(dat, significant_only = TRUE, ...) {
  if (significant_only) {
    significant <- abs(dat[["2.5%"]]) < 0
    dat_plot <- dat[significant, ]
  } else {
    dat_plot <- dat
  }
  if (nrow(dat_plot) > 0) {
    error_bar_plot(dat[["Mean"]], dat[["2.5%"]], dat[["97.5%"]], ...)
  } else {
    text_plot("Not significant")
  }
  box()
}

#' Simple text panel
#' 
#' @inheritParams graphics::text
text_plot <- function(x = 0, y = 0, labels, ...) {
  plot(0, 0, type = "n", axes = FALSE, xlab = NA, ylab = NA)
  text(x, y, labels = labels, ...)
}

#' Bar plot with error bars
#'
#' @inheritParams graphics::barplot
#' @param lower Error bar lower limit
#' @param upper Error bar upper limit
error_bar_plot <- function(height, lower, upper, ...) {
  centers <- barplot(height = height, ...)
  segments(centers, lower, centers, upper)
}

