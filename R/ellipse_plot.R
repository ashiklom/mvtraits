#' Retrieve covariance ellipse for single result
#'
#' @param mu List of mean coefficient statistics (see [fit_mvnorm()] and 
#' [fit_mvnorm_hier()])
#' @param Sigma List of variance-covariance matrix statistics (see 
#' [fit_mvnorm()] and [fit_mvnorm_hier()])
#' @inheritParams test_significant
#' @export
get_ellipse <- function(mu, Sigma, xvar, yvar) {
  xy <- c(xvar, yvar)
  mu_mean <- mu[["Mean"]][xy]
  Sigma_mean <- Sigma[["Mean"]][xy, xy]

  eps <- ellipse_axes(Sigma_mean, mu_mean)
  eps$significant <- test_significant(Sigma, xvar, yvar)
  eps
}

#' Retrieve all covariance ellipse data for plotting
#'
#' @param mu_global List of global mean statistics
#' @param Sigma_global List of global variance-covariance matrix statistics
#' @param mu_group List of group mean statistics
#' @param Sigma_group List of group variance-covariance matrix statistics
#' @param group_names Names of each group. Must match order in `mu_group` and 
#' `Sigma_group`.
#' @inheritParams get_ellipse
#' @export
prep_stickplot <- function(mu_global, Sigma_global, mu_group, Sigma_group,
                           xvar, yvar, group_names, global_name = "global") {
  xy <- c(xvar, yvar)

  eps_list <- Map(
    get_ellipse,
    mu = mu_group,  # apply over list
    Sigma = Sigma_group,  # apply over list
    xvar = xvar,
    yvar = yvar
  )
  eps_global <- get_ellipse(mu_global, Sigma_global, xvar, yvar)
  eps_dat <- dplyr::bind_rows(eps_global, eps_list)
  eps_dat$group <- c(global_name, group_names)
  eps_dat
}

#' Calculate x and y limits from ellipse data
#'
#' @param eps_dat Covariance ellipse data frame (see [prep_stickplot()])
#' @export
get_lims <- function(eps_dat) {
  xlim <- range(eps_dat$x1, eps_dat$x2)
  ylim <- range(eps_dat$y1, eps_dat$y2)
  list(xlim = xlim, ylim = ylim)
}

#' Draw ellipse covariance axes from ellipse data
#'
#' @inheritParams get_lims
#' @export
draw_sticks <- function(eps_dat) {
  cols <- colnames(eps_dat)
  stopifnot(c("center_x", "center_y", "x1", "x2", "y1", "y2") %in% cols)
  pch <- if ("pch" %in% cols) eps_dat$pch else 19
  col <- if ("col" %in% cols) eps_dat$col else "black"
  lwd <- if ("lwd" %in% cols) eps_dat$lwd else 1
  lty <- if ("lty" %in% cols) eps_dat$lty else "solid"
  cex <- if ("cex" %in% cols) eps_dat$cex else 1
  points(center_y ~ center_x, data = eps_dat,
         pch = pch, col = col, cex = cex)
  segments(eps_dat$x1, eps_dat$y1, eps_dat$x2, eps_dat$y2,
           col = col, lwd = lwd, lty = lty)
}

#' Stick plot for comparing group means and covariances
#'
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
default_stickplot <- function(mu_global, Sigma_global,
                              mu_group, Sigma_group,
                              xvar, yvar,
                              group_names, global_name = "global",
                              par_plot = list(), par_legend = list(),
                              col_global = "black",
                              col_group = RColorBrewer::brewer.pal(
                                length(group_names),
                                "Paired"
                              ),
                              xlab = xvar, ylab = yvar,
                              pch_g = 8, lty_g = 3, lwd_g = 3, cex_g = 2,
                              pch_s = 19, lty_s = 1, lwd_s = 2, cex_s = 1,
                              pch_n = 4, lty_n = 2, lwd_n = 0, cex_n = 1) {

  ngroup <- length(group_names)
  all_colors <- c(col_global, col_group)
  eps_dat <- prep_stickplot(
    mu_global, Sigma_global,
    mu_group, Sigma_group,
    xvar, yvar,
    group_names, global_name
  ) %>%
    dplyr::mutate(
      pch = dplyr::case_when(
        group == global_name ~ pch_g,
        significant ~ pch_s,
        TRUE ~ pch_n
      ),
      lty = dplyr::case_when(
        group == global_name ~ lty_g,
        significant ~ lty_s,
        TRUE ~ lty_n
      ),
      lwd = dplyr::case_when(
        group == global_name ~ lwd_g,
        significant ~ lwd_s,
        TRUE ~ lwd_n
      ),
      col = all_colors,
      cex = 1
    )
  lims <- get_lims(eps_dat)

  layout(matrix(c(1, 2), 1), widths = c(3, 1))
  par(par_plot)
  plot(0, 0, type = "n", xlim = lims$xlim, ylim = lims$ylim,
       xlab = "", ylab = "", axes = FALSE, ann = FALSE)
  draw_sticks(eps_dat)
  par(append(par_plot, list(new = TRUE)))
  plot(x = 1, y = 1, type = "n",
       xlim = 10 ^ lims$xlim, ylim = 10 ^ lims$ylim,
       log = "xy", xlab = xlab, ylab = ylab)

  par(append(par_legend, list(new = TRUE)))
  plot(x = 0, y = 0, type = "n", axes = FALSE, ann = FALSE)
  with(eps_dat, {
       legend("top", group, pch = pch, lty = lty,
              col = col, lwd = lwd, cex = cex)
       legend("bottom", c("significant", "not significant"),
              lty = c(lty_s, lty_n), pch = c(pch_s, pch_n))
       })
}
