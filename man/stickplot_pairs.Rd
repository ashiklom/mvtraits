% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ellipse_plot.R
\name{stickplot_pairs}
\alias{stickplot_pairs}
\title{Draw a pairs-like plot of covariance sticks}
\usage{
stickplot_pairs(mu_global_lower, Sigma_global_lower, mu_group_lower,
  Sigma_group_lower, vars_lower, group_names, vars_label = vars_lower,
  unlog_axes = FALSE, mu_global_upper = mu_global_lower,
  Sigma_global_upper = Sigma_global_lower,
  mu_group_upper = mu_group_lower,
  Sigma_group_upper = Sigma_group_lower, vars_upper = vars_lower,
  par_plot = list(), par_legplot = list(), par_legend = list(),
  par_label = list(), reorder_legend = FALSE,
  screen_split = rbind(c(0, 1, 0.2, 1), c(0, 1, 0, 0.2)), ...)
}
\arguments{
\item{mu_global_lower}{List of global means for lower triangle}

\item{Sigma_global_lower}{List of global variance-covariance matrices for lower triangle}

\item{mu_group_lower}{List of group means for lower triangle}

\item{Sigma_group_lower}{List of group variance-covariance matrices for lower triangle}

\item{vars_lower}{Names of variables for lower triangle}

\item{group_names}{Names of each group. Must match order in \code{mu_group} and
\code{Sigma_group}.}

\item{vars_label}{Labels for diagnoals}

\item{unlog_axes}{If TRUE, apply \code{10^x} transformation to axes}

\item{mu_global_upper}{List of global means for upper triangle}

\item{Sigma_global_upper}{List of global variance-covariance matrices for upper triangle}

\item{mu_group_upper}{List of group means for upper triangle}

\item{Sigma_group_upper}{List of group variance-covariance matrices for upper triangle}

\item{vars_upper}{Names of variables for upper triangle}

\item{par_plot}{List of graphical parameters for plot}

\item{par_legplot}{List of graphical parameters for legend plot box}

\item{par_legend}{List of options passed directly to \code{legend} function}

\item{par_label}{List of graphical parameters for diagonal labels}

\item{reorder_legend}{If \code{TRUE}, put last drawn layer first in legend.}

\item{screen_split}{\code{figs} argument for \code{\link[graphics:split_screen]{graphics::split_screen()}} for plot-legend split.}

\item{...}{Additional parameters to \code{\link[=prep_stickplot]{prep_stickplot()}}}
}
\description{
Draw a pairs-like plot of covariance sticks
}
