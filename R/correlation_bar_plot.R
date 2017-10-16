#' Draw matrix of correlation bar plots
#'
#' @param dat_lower Data for lower triangle
#' @param dat_upper Data for upper triangle
#' @param param_order Character vector of parameter values indicating order
#' @param ... Additional arguments to [single_bar_plot]
bar_plot_matrix <- function(dat_lower, dat_upper,
                            param_order = NULL, significant_only = TRUE,
                            diag_cex = 2, ...) {
  if (is.null(param_order)) {
    param_order <- union(dat_lower[["xvar"]], dat_lower[["yvar"]])
  }
  xind_lower <- xind_upper <- match(dat_lower[["xvar"]], param_order)
  yind_lower <- yind_upper <- match(dat_lower[["yvar"]], param_order)
  dat_lower_r <- dat_lower[order(yind_lower, xind_lower), ]
  dat_upper_r <- dat_upper[order(yind_upper, xind_upper), ]

  nparam <- length(param_order)
  npseq <- seq_len(0.5 * (nparam * (nparam - 1)))
  npseq_max <- max(npseq)
  lt <- ut <- matrix(0, nparam, nparam)
  lt[lower.tri(lt)] <- npseq
  ut[lower.tri(ut)] <- npseq_max + npseq
  lay_mat <- lt + t(ut)
  diag(lay_mat) <- 2 * npseq_max + seq_len(nparam)
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
  lapply(param_order, text_plot, cex = diag_cex)
}

#' Single bar plot panel
#'
#' @param dat Data for panel
#' @param significant_only If `TRUE`, hide columns that are not significantly
#' difference from zero.
single_bar_plot <- function(dat, significant_only = TRUE, ...) {
  if (significant_only) {
    significant <- sign(dat[["2.5%"]]) == sign(dat[["97.5%"]])
    dat_plot <- dat[significant, ]
  } else {
    dat_plot <- dat
  }
  if (nrow(dat_plot) > 0) {
    error_bar_plot(dat_plot[["Mean"]], dat_plot[["2.5%"]], dat_plot[["97.5%"]], ...)
  } else {
    text_plot("Not significant")
  }
  box()
}

#' Simple text panel
#'
#' @inheritParams graphics::text
text_plot <- function(labels, x = 0, y = 0, ...) {
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

