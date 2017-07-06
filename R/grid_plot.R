#' Column ggplot of summary statistics
#' 
#' @param summary_data `tibble` containing summary statistics (e.g. output of `summary_df`)
#' @param varname Character string of variable to plot (e.g. "mu", "Sigma", "Corr")
#' @param aes_list List of parameters passed to `aes_string`
#' @param col_list List of parameters passed to `geom_col`
#' @param errorbar_list List of parameters passed to `geom_errorbar`
#' @param facet_fun Function for facetting. Default is `facet_wrap`. Use `NULL` to turn off facetting.
#' @param facet_list List of parameters passed to `facet_wrap`
#' @param theme_list List of parameters passed to `theme`.
#' @export
summary_barplot <- function(summary_data, varname, 
                         aes_list = list(),
                         errorbar_list = list(),
                         col_list = list(),
                         facet_fun = ggplot2::facet_wrap,
                         facet_list = list(),
                         theme_list = list()) {
    aes_list <- modifyList(list(x = 'pft', y = 'Mean', 
                                ymin = '`2.5%`', ymax = '`97.5%`', 
                                fill = 'pft'),
                           aes_list)
    col_list = modifyList(list(), col_list)
    errorbar_list = modifyList(list(width = 0.2), errorbar_list)
    facet_list = modifyList(list(facets = formula(~xparam + yparam), drop = FALSE), facet_list)
    theme_list <- modifyList(list(axis.text.x = ggplot2::element_blank(),
                                  axis.ticks.x = ggplot2::element_blank()), 
                             theme_list)
    dat_plot <- dplyr::filter(summary_data, variable == varname)
    plt <- ggplot2::ggplot(dat_plot) + 
        do.call(ggplot2::aes_string, aes_list) +
        do.call(ggplot2::geom_col, col_list) + 
        do.call(ggplot2::geom_errorbar, errorbar_list)
    if (!is.null(facet_fun)) {
        plt <- plt + do.call(facet_fun, facet_list)
    }
    plt <- plt + do.call(ggplot2::theme, theme_list)
    return(plt)
}
