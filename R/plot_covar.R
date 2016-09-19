#' @export
plot_covar <- function(sumtab, ...){
    library(corrplot)
    cor_global <- tab2mat(sumtab, corr = TRUE)
    plt <- corrplot.mixed(cor_global, lower = "ellipse", upper = "number",
                          mar=c(0,0,2,0), ...)
    return(plt)
}
