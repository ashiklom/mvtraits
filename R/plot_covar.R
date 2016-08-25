plot_covar <- function(fname, ...){
    library(corrplot)
    cov.global <- readRDS(fname)$Sigma_trait
    cor.global <- cov.global %>%
        covToCor.global %>%
        colMeans %>%
        matrix(sqrt(length(.)))
    plt <- corrplot.mixed(cor.global, lower = "ellipse", upper = "number",
                          mar=c(0,0,2,0), ...)
    return(plt)
}
