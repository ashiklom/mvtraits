#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom rlang UQ UQS

#' @export
summary2vec <- function(summary_table, val_col = "Mean", ...) {
  filterq <- rlang::quos(...)
  tab_sub <- dplyr::filter(summary_table, UQS(filterq))
  vals <- tab_sub[[val_col]]
  names(vals) <- tab_sub[["index"]]
  vals
}

#' @export
summary2mat <- function(summary_table, val_col = "Mean", ...) {
  filterq <- rlang::quos(...)
  tab_sub <- summary_table %>%
    dplyr::filter(UQS(filterq)) %>%
    tidyr::separate(index, c("var1", "var2"), sep = "\\.\\.") %>%
    dplyr::select(var1, var2, dplyr::one_of(val_col))
  if (is.factor(tab_sub$var1)) {
    params <- levels(tab_sub$var1)
  } else {
    params <- union(unique(tab_sub$var1), unique(tab_sub$var2))
  }
  dfcol2mat(
    column = tab_sub[[val_col]],
    xvar = tab_sub[["var1"]],
    yvar = tab_sub[["var2"]],
    params = params
  )
}
