#' Prepare covariance ellipse for plotting
#'
#' @param mu List of mean coefficient statistics (see [fit_mvnorm()] and 
#' [fit_mvnorm_hier()])
#' @param Sigma List of variance-covariance matrix statistics (see 
#' [fit_mvnorm()] and [fit_mvnorm_hier()])
#' @inheritParams test_significant
#' @param lty_s Line type for significant correlations
#' @param lty_n Line type for insignificant correlations
#' @param lwd_s Line width for significant correlations
#' @param lwd_n Line width for insignificant correlations
#' @param pch_s Point shape for significant correlations
#' @param pch_n Point shape for insignificant correlations
#' @param hide_insignificant If TRUE, set coordinates for insignificant 
#' correlations to `NA`
prep_ellipse <- function(mu, Sigma,
                          xvar, yvar,
                          lty_s = "solid", lty_n = "dashed",
                          lwd_s = 2.0, lwd_n = 0.5,
                          pch_s = 19, pch_n = 4,
                          hide_insignificant = FALSE) {
  xy <- c(xvar, yvar)
  mu_mean <- mu[["Mean"]][xy]
  Sigma_mean <- Sigma[["Mean"]][xy, xy]

  is_significant_global <- test_significant(Sigma, xvar, yvar)
  eps <- ellipse_axes(Sigma_mean, mu_mean)
  if (is_significant_global) {
    eps$lty <- lty_s
    eps$lwd <- lwd_s
    eps$pch <- pch_s
  } else {
    if (hide_insignificant) {
      eps$x1 <- eps$x2 <- eps$y1 <- eps$y2 <- NA
    }
    eps$lty <- lty_n
    eps$lwd <- lwd_n
    eps$pch <- pch_n
  }
  eps
}

#' Stick plot for comparing group means and covariances
#'
#' @param mu_global List of global mean statistics
#' @param Sigma_global List of global variance-covariance matrix statistics
#' @param mu_group List of group mean statistics
#' @param Sigma_group List of group variance-covariance matrix statistics
#' @param group_names Names of each group. Must match order in `mu_group` and 
#' `Sigma_group`.
#' @param par_plot List of graphical parameters for plot
#' @param par_legend List of graphical parameters for legend
#' @param cols Colors for group values
#' @param col_g Color for global values
#' @param xlab Label for x variable
#' @param ylab Label for y variable
#' @param hide_insignificant_global If TRUE, don't draw global correlation if 
#' it isn't statistically significant
#' @param lty_g Line type for global correlation
#' @param lwd_g Line width for global correlation
#' @param pch_g Point shape for global correlation
#' @param cex_g Point size for global correlation
#' @inheritParams prep_ellipse
#' @export
stick_plot <- function(mu_global, Sigma_global,
                       mu_group, Sigma_group,
                       xvar, yvar, group_names,
                       par_plot = par(), par_legend = par(),
                       cols = RColorBrewer::brewer.pal(length(group_names), 
                                                       "Paired"),
                       xlab = xvar, ylab = yvar,
                       hide_insignificant_global = FALSE,
                       col_g = "black",
                       pch_g = 8, lty_g = "dotted", lwd_g = 3, cex_g = 2,
                       hide_insignificant = TRUE,
                       lty_s = "solid", lty_n = "dashed",
                       lwd_s = 2.0, lwd_n = 0.5,
                       pch_s = 19, pch_n = 4) {
  xy <- c(xvar, yvar)
  ngroup <- length(group_names)

  eps_list <- Map(
    prep_ellipse,
    mu = mu_group,  # apply over list
    Sigma = Sigma_group,  # apply over list
    xvar = xvar,
    yvar = yvar,
    hide_insignificant = hide_insignificant,
    lty_s = lty_s, lty_n = lty_n,
    lwd_s = lwd_s, lwd_n = lwd_n,
    pch_s = pch_s, pch_n = pch_n
  )
  eps_dat <- dplyr::bind_rows(eps_list)
  eps_dat$group <- group_names

  eps_global <- ellipse_axes(
    Sigma_global$Mean[xy, xy],
    mu_global$Mean[xy]
  )
  if (hide_insignificant_global) {
    is_significant_global <- test_significant(Sigma_global, xvar, yvar)
  } else {
    is_significant_global <- TRUE
  }

  xlim <- with(eps_dat, c(min(x1, x2, na.rm = TRUE), max(x1, x2, na.rm = TRUE)))
  ylim <- with(eps_dat, c(min(y1, y2, na.rm = TRUE), max(y1, y2, na.rm = TRUE)))

  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(par_plot)
  plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
       xlab = "", ylab = "", axes = FALSE)

  ## Group points
  with(eps_dat, points(center_x, center_y, pch = pch, col = cols))
  with(eps_dat, segments(x1, y1, x2, y2, col = cols, lwd = lwd, lty = lty))

  ## Global correlation
  with(eps_global, points(center_x, center_y,
                          pch = pch_g, col = col_g, cex = cex_g))
  if (is_significant_global) {
    with(eps_global, segments(x1, y1, x2, y2,
                              col = col_g, lwd = lwd_g, lty = lty_g))
  }

  labels_l <- c("global", group_names)
  pch_l <- c(pch_g, rep(pch_s, ngroup))
  lty_l <- c(lty_g, rep(lty_s, ngroup))
  col_l <- c(col_g, cols)

  ## Add log10 axis
  par(append(par_plot, list(new = TRUE)))
  plot(1, 1, type = "n", xlim = 10 ^ xlim, ylim = 10 ^ ylim,
       log = "xy", xlab = xlab, ylab = ylab)

  par(par_legend)
  plot(0:3, 0:3, type = "n", ann = FALSE, axes = FALSE)
  legend("topleft", labels_l, pch = pch_l, lty = lty_l, col = col_l, 
         bty = "n")

  if (hide_insignificant) {
    lty_ln <- NA
  } else {
    lty_ln <- lty_n
  }

  legend("bottomleft", c("significant", "not significant"),
         lty = c(lty_s, lty_ln), pch = c(pch_s, pch_n), bty = "n")
}
