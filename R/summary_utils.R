#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom rlang UQ

summary2vec_multi <- function(summary_table, variable = "mu") {
  tab_sub <- dplyr::filter(summary_table, variable == UQ(variable))
  mean_vals <- tab_sub[["Mean"]]
  names(mean_vals) <- tab_sub[["index"]]
  mean_vals
}

summary2mat_multi <- function(summary_table, variable = 'Sigma') {
  tab_sub <- summary_table %>% 
    dplyr::filter(variable == UQ(variable)) %>% 
    tidyr::separate(index, c("var1", "var2"), sep = "\\.\\.") %>% 
    dplyr::select(var1, var2, Mean)
  if (is.factor(tab_sub$var1)) {
    params <- levels(tab_sub$var1)
  } else {
    params <- union(unique(tab_sub$var1), unique(tab_sub$var2))
  }
  with(tab_sub, dfcol2mat(column = Mean, xvar = var1, yvar = var2, params = params))
}

summary2mu_global <- function(summary_table) {
  tab_sub <- dplyr::filter(summary_table, variable == 'mu')
  mean_vals <- tab_sub[["Mean"]]
  names(mean_vals) <- tab_sub[["index"]]
  mean_vals
}

summary2mu_group <- function(summary_table) {
  tab_sub <- dplyr::filter(summary_table, variable == 'mu')
  mean_vals <- tab_sub[["Mean"]]
  names(mean_vals) <- tab_sub[["index"]]
  mean_vals
}
