#' @importFrom magrittr `%>%`
#' @export
summary_df <- function(results_mcmc_list, group = NULL) {
    stopifnot(coda::is.mcmc.list(results_mcmc_list))
    rawsum <- summary(results_mcmc_list)
    procsum <- purrr::map(c('statistics', 'quantiles'), ~tibble_ize(rawsum[[.]], 'rowname'))
    dat_merged <- Reduce(dplyr::full_join, procsum)
    if (isTRUE(group)) {
        result <- dat_merged %>% 
            tidyr::separate(rowname, c('variable', 'group', 'index'), 
                            sep = varsep_esc, remove = FALSE, extra = 'merge')
    } else {
        result <- dat_merged %>% 
            tidyr::separate(rowname, c('variable', 'index'), 
                            sep = varsep_esc, remove = FALSE, extra = 'merge')
        if (!is.null(group)) {
            result <- dplyr::mutate(result, group = group)
        }
    }
    return(result)
}

tibble_ize <- function(mat, var) {
    as.data.frame(mat) %>% 
        tibble::rownames_to_column(var = var) %>% 
        tibble::as_tibble()
}
