#' @export
plot_covar <- function(sumtab, ...){
    library(corrplot)
    dots <- list(...)
    tab2mat_args <- c("colorder")
    colorder <- dots[["colorder"]]
    dots_plot <- dots[!names(dots) %in% tab2mat_args]
    cor_global <- tab2mat(sumtab, corr = TRUE, colorder = colorder)
    plt <- do.call(corrplot.mixed, c(list(corr = cor_global,
                                          lower = "ellipse", 
                                          upper = "number",
                                          mar = c(0,0,2,0)),
                                     dots_plot))
    return(plt)
}
