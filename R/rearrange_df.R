#' Rearrange result data frame into correct parameter order
#'
#' @param dat_full Summary tibble
#' @param params Vector of params, in correct order
#' @export
rearrange_df <- function(dat_full, params) {
    nparams <- length(params)
    dat_bignest <- tidyr::nest(dat_full, -variable)
    dat_bigproc <- dplyr::mutate(dat_bignest,
                                 data = purrr::map2(data, variable, ~procdat(.x, params, .y)))
    dat_unnest <- tidyr::unnest(dat_bigproc)
    return(dat_unnest)
}

procdat <- function(dat, params, variable) {
    if (variable == 'mu') {
        out <- proc_vec(dat, params)
    } else if (variable %in% c('Sigma', 'Corr')) {
        out2 <- proc_mat(dat, params)
        if (variable == 'Corr') {
            out <- filter(out2, as.character(xparam) != as.character(yparam))
        } else {
            out <- out2
        }
    }
    return(out)
}

proc_vec <- function(dat, params) {
    dat %>%
        dplyr::mutate(param = factor(index, levels = params)) %>%
        dplyr::select(-index)
}

proc_mat <- function(dat, params) {
    oldnames <- names(dat)
    newnames <- make.names(oldnames)
    colnames(dat) <- newnames
    value_cols <- oldnames[!oldnames %in% c('group', 'index')]
    dat_nested <- tidyr::nest(dat, -group)
    out_nested <- dplyr::mutate(dat_nested,
                                data = purrr:::map(data, reorder_df,
                                                   params = params, value_cols = value_cols))
    out <- tidyr::unnest(out_nested)
    return(out)
}

reorder_df <- function(dat, params, value_cols) {
    dat_split <- tidyr::separate(dat, index, c('xvar', 'yvar'), sep = varsep_esc)
    out <- dat_split %>%
        dplyr::select_if(is.double) %>%
        lapply(dfcol2df, xvar = dat_split[['xvar']], yvar = dat_split[['yvar']], params) %>%
        Reduce(f = function(x, y) dplyr::full_join(x, y, by = c('xparam', 'yparam')), x = .)
    colnames(out)[grep('value', colnames(out))] <- value_cols
    return(out)
}

dfcol2df <- function(column, xvar, yvar, params) {
    mat <- dfcol2mat(column, xvar, yvar, params)
    vec <- flatten_matrix(mat)
    new_vars <- strsplit(names(vec), split = varsep_esc)
    xparam <- sapply(new_vars, '[[', 1)# %>% factor(levels = params[-1])
    yparam <- sapply(new_vars, '[[', 2)# %>% factor(levels = params[-length(params)])
    out <- tibble::tibble(xparam = xparam, yparam = yparam, value = vec)
    return(out)
}

#' @export
dfcol2mat <- function(column, xvar, yvar, params) {
    m <- length(params)
    n <- length(column)
    mat <- diag(1, m)
    dimnames(mat) <- list(params, params)
    for (i in seq_len(n)) {
        xx <- xvar[i]
        yy <- yvar[i]
        mat[xx, yy] <- column[i]
        mat[yy, xx] <- column[i]
    }
    return(mat)
}

