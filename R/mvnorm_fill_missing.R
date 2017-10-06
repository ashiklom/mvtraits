#' @export
setup_missing <- function(dat) {
    # Use -1 throughout to make indices work with C code
    all_missing <- which(apply(dat, 2, function(x) all(is.na(x)))) - 1
    missing_vec <- apply(dat, 1, function(x) which(is.na(x)) - 1)
    missing_vec_unique <- unique(missing_vec)
    missing_pattern <- sapply(missing_vec, paste, collapse = '_')
    missing_pattern_unique <- unique(missing_pattern)
    pattern_inds <- sapply(missing_pattern_unique, function(x) which(missing_pattern == x) - 1)
    present_vec <- apply(dat, 1, function(x) which(!is.na(x)) - 1)
    present_vec_unique <- unique(present_vec)
    n_pattern <- length(missing_pattern_unique)
    out <- list(npatt = n_pattern, mlist = missing_vec_unique, plist = present_vec_unique,
                indlist = pattern_inds, all_missing = all_missing)
    return(out)
}

#' @export
mvnorm_fill_missing <- function(dat, mu, Sigma, setup = NULL) {
    # Split data into missing pieces
    if (is.null(setup)) {
        setup <- setup_missing(dat)
    }
    dat_filled <- c_mvnorm_fill_missing(dat, mu, Sigma, setup)
    return(dat_filled)
}
