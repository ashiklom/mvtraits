#' @export
pairs_density <- function(uni.mus, multi.mus, hier.mus, 
                          obs.means, obs.global=NA, ...){

    nsamp <- min(nrow(uni.mus), nrow(multi.mus), nrow(hier.mus), 5000)

    uni <- 1:nsamp
    multi <- uni + nsamp
    hier <- multi + nsamp
    breaks.dens <- list(uni, multi, hier)

    obs_mean <- nsamp*3 + 1
    global_mean <- obs_mean + 1
    breaks.line <- list(obs_mean, global_mean)

    breaks <- list("density" = breaks.dens, "line" = breaks.line)

    dens.pars <- list(uni = list(col = "blue", lwd=2),
                      multi = list(col = "green4", lwd=2),
                      hier = list(col = "red", lwd=2))

    line.pars <- list(obs_means = list(lty = "dashed",
                                       col = "black", 
                                       lwd = 2),
                      global_means = list(lty = "dotted", 
                                          col="hotpink4", 
                                          lwd = 2))

    mus.list <- list(uni.mus[sample.int(nsamp),],
                     multi.mus[sample.int(nsamp),],
                     hier.mus[sample.int(nsamp),],
                     obs.means, 
                     obs.global)
    plot.mus <- do.call(rbind, mus.list)

    pairs(plot.mus, panel = dens2d.panel, densfunc = ellipse2d,
          breaks = breaks,
          dens.pars = dens.pars, line.pars = line.pars, ...)
}

#' @export
dens2d <- function(x,y,...) {
    z <- MASS::kde2d(x, y)
    levs <- c(0.8)
    contour(z, drawlabels=FALSE, levels=levs, add=TRUE, ...)
}

#' @export
ellipse2d <- function(x, y, ...){
    ellipse::dataEllipse(x, y, levels = 0.975, add=TRUE, 
                         plot.points = FALSE, center.pch=FALSE, ...)
}

#' @export
dens2d.panel <- function(x, y, densfunc, breaks, dens.pars, line.pars, ...){
    breaks.dens <- breaks[["density"]]
    breaks.line <- breaks[["line"]]

    for(i in seq_along(breaks.dens)){
        ind <- breaks.dens[[i]]
        plot_args <- c(list(x = x[ind], y = y[ind]), dens.pars[[i]])
        do.call(densfunc, plot_args)
    }

    for(i in seq_along(breaks.line)){
        ind <- breaks.line[[i]]
        plot_args <- c(list(v = x[ind], h = y[ind]), line.pars[[i]])
        do.call(abline, plot_args)
    }
}
