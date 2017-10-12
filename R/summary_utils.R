#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom rlang UQ UQS

summary2vec <- function(summary_table, ...) {
  filterq <- rlang::quos(...)
  tab_sub <- dplyr::filter(summary_table, UQS(filterq))
  mean_vals <- tab_sub[["Mean"]]
  names(mean_vals) <- tab_sub[["index"]]
  mean_vals
}

summary2mat <- function(summary_table, ...) {
  filterq <- rlang::quos(...)
  tab_sub <- summary_table %>% 
    dplyr::filter(UQS(filterq)) %>% 
    tidyr::separate(index, c("var1", "var2"), sep = "\\.\\.") %>% 
    dplyr::select(var1, var2, Mean)
  if (is.factor(tab_sub$var1)) {
    params <- levels(tab_sub$var1)
  } else {
    params <- union(unique(tab_sub$var1), unique(tab_sub$var2))
  }
  with(tab_sub, dfcol2mat(column = Mean, xvar = var1, yvar = var2, params = params))
}
