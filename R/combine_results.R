combine_results <- function(prev_results, curr_results) {
    stopifnot(length(prev_results) == length(curr_results),
              names(prev_results[[1]]) == names(curr_results)[[1]])
    nchain <- length(prev_results)
    vars <- names(prev_results[[1]])
    out <- list()
    for (i in seq_len(nchain)) {
        out[[i]] <- list()
        for (v in vars) {
            out[[i]][[v]] <- abind::abind(prev_results[[i]][[v]], curr_results[[i]][[v]], along = 1)
        }
    }
    return(out)
}
