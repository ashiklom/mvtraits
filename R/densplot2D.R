dens2d <- function(x,y,...) {
    z <- kde2d(x, y)
    levs <- c(0.8)
    contour(z, drawlabels=FALSE, levels=levs, add=TRUE, ...)
}

ellipse2d <- function(x, y, ...){
    dataEllipse(x, y, levels = 0.975, add=TRUE, 
                plot.points = FALSE, center.pch=FALSE, ...)
}

dens2d.panel <- function(x, y, densfunc, breaks, dens.pars, line.pars){
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
