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
#' @param global_name Name of global group. Default = "global".
#' @inheritParams get_ellipse
#' @param col_global Color for global correlation
#' @param col_group List of colors for group correlations
#' @param lty_g Line type for global correlation
#' @param lwd_g Line width for global correlation
#' @param pch_g Point shape for global correlation
#' @param cex_g Point size for global correlation
#' @param lty_s Line type for significant group correlations
#' @param lwd_s Line width for significant group correlations
#' @param pch_s Point shape for significant group correlations
#' @param cex_s Point size for significant group correlations
#' @param lty_n Line type for insignificant group correlations
#' @param lwd_n Line width for insignificant group correlations
#' @param pch_n Point shape for insignificant group correlations
#' @param cex_n Point size for insignificant group correlations
#' @export
prep_stickplot <- function(mu_global, Sigma_global, mu_group, Sigma_group,
                           xvar, yvar, group_names, global_name = "global",
                           col_global = "black",
                           col_group = RColorBrewer::brewer.pal(
                             length(group_names), "Paired"
                           ),
                           pch_g = 8, lty_g = 3, lwd_g = 3, cex_g = 2,
                           pch_s = 19, lty_s = 1, lwd_s = 2, cex_s = 1,
                           pch_n = 4, lty_n = 0, lwd_n = 1, cex_n = 1
                           ) {
  xy <- c(xvar, yvar)

  eps_list <- Map(
    get_ellipse,
    mu = mu_group,  # apply over list
    Sigma = Sigma_group,  # apply over list
    xvar = xvar,
    yvar = yvar
  )

  set_pars <- function(dat, g, s, n, global_name = "global", global_pri = FALSE) {
    if (global_pri) {
      dplyr::case_when(
        dat$group == global_name ~ g,
        !dat$significant ~ n,
        dat$significant ~ s
      )
    } else {
      dplyr::case_when(
        !dat$significant ~ n,
        dat$group == global_name ~ g,
        dat$significant ~ s
      )
    }
  }

  eps_global <- get_ellipse(mu_global, Sigma_global, xvar, yvar)
  eps_dat <- dplyr::bind_rows(eps_list, eps_global)
  eps_dat$group <- c(group_names, global_name)
  all_colors <- c(col_group, col_global)
  eps_dat %>%
    dplyr::mutate(
      pch = set_pars(.data, pch_g, pch_s, pch_n, global_name, TRUE),
      lty = set_pars(.data, lty_g, lty_s, lty_n, global_name),
      lwd = set_pars(.data, lwd_g, lwd_s, lwd_n, global_name),
      cex = set_pars(.data, cex_g, cex_s, cex_n, global_name),
      col = all_colors,
      # store defaults for legend
      pch_l = dplyr::if_else(group == global_name, pch_g, pch_s),
      lty_l = dplyr::if_else(group == global_name, lty_g, lty_s),
      lwd_l = dplyr::if_else(group == global_name, lwd_g, lwd_s)
    )
}

#' Calculate x and y limits from ellipse data
#'
#' @param eps_dat Covariance ellipse data frame (see [prep_stickplot()])
#' @export
get_lims <- function(eps_dat) {
  eps_lim <- dplyr::mutate_at(
    eps_dat,
    c("x1", "x2", "y1", "y2"),
    dplyr::funs(dplyr::if_else(significant, ., NA_real_))
  )
  xlim <- range(eps_lim$center_x, eps_lim$x1, eps_lim$x2, na.rm = TRUE)
  ylim <- range(eps_lim$center_y, eps_lim$y1, eps_lim$y2, na.rm = TRUE)
  list(xlim = xlim, ylim = ylim)
}

#' Draw ellipse covariance axes from ellipse data
#'
#' @inheritParams get_lims
#' @export
draw_sticks <- function(eps_dat) {
  cols <- colnames(eps_dat)
  stopifnot(c("center_x", "center_y", "x1", "x2", "y1", "y2") %in% cols)
  for (i in seq_len(nrow(eps_dat))) {
    d <- eps_dat[i,]
    pch <- if ("pch" %in% cols) d$pch else 19
    col <- if ("col" %in% cols) d$col else "black"
    lwd <- if ("lwd" %in% cols) d$lwd else 1
    lty <- if ("lty" %in% cols) d$lty else "solid"
    cex <- if ("cex" %in% cols) d$cex else 1
    segments(d$x1, d$y1, d$x2, d$y2,
             col = col, lwd = lwd, lty = lty)
    points(d$center_x, d$center_y, pch = pch, col = col, cex = cex)
  }
}

#' Stick plot for comparing group means and covariances
#'
#' @inheritParams prep_stickplot
#' @param par_plot List of graphical parameters for plot
#' @param par_legend List of graphical parameters for legend
#' @param xlab Label for x variable
#' @param ylab Label for y variable
#' @param unlog_axes If TRUE, apply `10^x` transformation to axes
#' @param ... Additional parameters to [prep_stickplot()]
#' @export
default_stickplot <- function(mu_global, Sigma_global, mu_group, Sigma_group,
                              xvar, yvar, group_names,
                              par_plot = list(), par_legend = list(),
                              xlab = xvar, ylab = yvar,
                              unlog_axes = FALSE,
                              ...) {

  eps_dat <- prep_stickplot(
    mu_global, Sigma_global,
    mu_group, Sigma_group,
    xvar, yvar, group_names,
    ...
  )
  lims <- get_lims(eps_dat)

  f1 <- c(0, 0.75, 0, 1)
  f2 <- c(0.75, 1, 0, 1)
  par(modifyList(list(fig = f1), par_plot))
  plot(0, 0, type = "n", xlim = lims$xlim, ylim = lims$ylim,
       xlab = "", ylab = "", axes = FALSE, ann = FALSE)
  draw_sticks(eps_dat)
  if (unlog_axes) {
    par(modifyList(list(new = TRUE, fig = f1), par_plot))
    plot(x = 1, y = 1, type = "n",
         xlim = 10 ^ lims$xlim, ylim = 10 ^ lims$ylim,
         log = "xy", xlab = xlab, ylab = ylab)
  } else {
    box(which = "plot", lty = "solid", col = "black")
    axis(side = 1)
    mtext(xlab, side = 1, line = 3)
    axis(side = 2)
    mtext(ylab, side = 2, line = 3)
  }

  defpar_legend <- list(
    new = TRUE,
    fig = f2,
    mar = c(5.1, 0.1, 4.1, 0.1)
  )
  par(modifyList(defpar_legend, par_legend))
  plot(0, 0, type = "n", axes = FALSE, ann = FALSE)
  with(eps_dat, {
       legend("topleft", group, 
              pch = pch, lty = lty, col = col, lwd = lwd,
              bty = "n")
       })
}

#' Draw a pairs-like plot of covariance sticks
#'
#' @param mu_global_lower List of global means for lower triangle
#' @param Sigma_global_lower List of global variance-covariance matrices for lower triangle
#' @param mu_group_lower List of group means for lower triangle
#' @param Sigma_group_lower List of group variance-covariance matrices for lower triangle
#' @param vars_lower Names of variables for lower triangle
#' @inheritParams default_stickplot
#' @param vars_label Labels for diagnoals
#' @param mu_global_upper List of global means for upper triangle
#' @param Sigma_global_upper List of global variance-covariance matrices for upper triangle
#' @param mu_group_upper List of group means for upper triangle
#' @param Sigma_group_upper List of group variance-covariance matrices for upper triangle
#' @param vars_upper Names of variables for upper triangle
#' @param par_label List of graphical parameters for diagonal labels
#' @param par_legplot List of graphical parameters for legend plot box
#' @param par_legend List of options passed directly to `legend` function
#' @param reorder_legend If `TRUE`, put last drawn layer first in legend.
#' @param screen_split `figs` argument for [graphics::split_screen()] for plot-legend split.
#' @param ... Additional parameters to [prep_stickplot()]
#' @export
stickplot_pairs <- function(mu_global_lower, Sigma_global_lower,
                            mu_group_lower, Sigma_group_lower,
                            vars_lower, group_names,
                            vars_label = vars_lower, unlog_axes = FALSE,
                            mu_global_upper = mu_global_lower,
                            Sigma_global_upper = Sigma_global_lower,
                            mu_group_upper = mu_group_lower,
                            Sigma_group_upper = Sigma_group_lower,
                            vars_upper = vars_lower,
                            par_plot = list(), 
                            par_legplot = list(),
                            par_legend = list(),
                            par_label = list(),
                            reorder_legend = FALSE,
                            screen_split = rbind(c(0, 1, 0.2, 1), c(0, 1, 0, 0.2)),
                            ...) {
  stopifnot(length(vars_lower) == length(vars_upper))
  nv <- length(vars_lower)
  main_screens <- split.screen(screen_split)
  par(modifyList(list(mar = rep(0.1, 4), oma = c(0, 3, 3, 3)), par_plot))
  screens <- split.screen(c(nv, nv), screen = main_screens[1])
  for (i in seq_len(nv)) {
    for (j in seq_len(nv)) {
      k <- j + nv * (i - 1)
      screen(screens[k])
      if (i == j) {
        # Label
        plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
        box()
        text_default <- list(
          x = 0.5, y = 0.5, labels = vars_label[i]
        )
        do.call(text, modifyList(text_default, par_label))
      } else {
        if (j > i) {
          # Upper triangle
          mu_global <- mu_global_upper
          Sigma_global <- Sigma_global_upper
          mu_group <- mu_group_upper
          Sigma_group <- Sigma_group_upper
          xvar <- vars_upper[j]
          yvar <- vars_upper[i]
        } else {
          # Lower triangle
          mu_global <- mu_global_lower
          Sigma_global <- Sigma_global_lower
          mu_group <- mu_group_lower
          Sigma_group <- Sigma_group_lower
          xvar <- vars_lower[j]
          yvar <- vars_lower[i]
        }
        eps_dat <- prep_stickplot(
          mu_global, Sigma_global, mu_group, Sigma_group,
          xvar, yvar, group_names, ...
        )
        lims <- get_lims(eps_dat)
        plot(0, 0, type = "n", xlim = lims$xlim, ylim = lims$ylim,
             xlab = "", ylab = "", axes = FALSE, ann = FALSE)
        draw_sticks(eps_dat)
        if (unlog_axes) {
          screen(screens[k])
          plot(x = 1, y = 1, type = "n",
               xlim = 10 ^ lims$xlim, ylim = 10 ^ lims$ylim,
               log = "xy", xlab = "", ylab = "", axes = FALSE, ann = FALSE)
        }
        box(which = "plot", lty = "solid", col = "black")
        if (i == 1) {
          axis(side = 3)
        }
        if (j == 1) {
          axis(side = 2)
        }
        if (i == nv) {
          axis(side = 1)
        }
        if (j == nv) {
          axis(side = 4)
        }
      }
    }
  }
  ## Draw legend
  screen(main_screens[2])
  par(par_legplot)
  plot(0, 0, type = "n", axes = FALSE, ann = FALSE)
  if (reorder_legend) {
    s <- c(nrow(eps_dat), seq(1, nrow(eps_dat) - 1))
    eps_leg <- eps_dat[s,]
  } else {
    eps_leg <- eps_dat
  }
  leg_default <- with(eps_leg, list(
    x = "center",
    legend = group,
    col = col,
    pch = pch_l,
    lty = lty_l,
    lwd = lwd_l,
    bty = "n",
    ncol = 4
    )
  )
  do.call(legend, modifyList(leg_default, par_legend))
  close.screen(all.screens = TRUE)
}
