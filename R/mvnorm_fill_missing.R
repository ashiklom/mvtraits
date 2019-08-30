#' Fill missing multivariate data based on known mean vector and covariance
#' matrix
#'
#' @param dat Data matrix (observations in rows, variables in columns) with
#'     missing data as NA
#' @param mu Mean vector. Must be length `ncol(dat)`
#' @param Sigma Covariance matrix. Must be dimensions `ncol(dat) x ncol(dat)`
#' @param setup Information about missingness pattern, as returned by
#'     [setup_missing()]. If missing, it will be calculated. When running this
#'     function many times on data with the same missingness pattern, run
#'     [setup_missing()] once first to make this calculation more efficient.
#' @return One realization of `dat` with missing values imputed
#' @export
mvnorm_fill_missing <- function(dat, mu, Sigma, setup = NULL) {
  # If everything in data is missing, just draw randomly from mu and sigma
  # (with a warning)
  if (all(is.na(dat))) {
    warning("Everything in `dat` is missing. Filling with random draws.")
    out <- mvtraits::random_mvnorm(nrow(dat), mu, Sigma)
    return(out)
  }

  # Split data into missing pieces
  if (is.null(setup)) {
    setup <- setup_missing(dat)
  }
  dat_filled <- c_mvnorm_fill_missing(dat, mu, Sigma, setup)
  return(dat_filled)
}

#' Analyze missingness pattern in data
#'
#' Used as an input for [mvnorm_fill_missing()].
#'
#' @inheritParams mvnorm_fill_missing
#' @return List describing the missingness pattern in `dat` and other relevant
#'     indices for sorting, etc.
#' @export
setup_missing <- function(dat) {
  # Use -1 throughout to make indices work with C code

  all_missing <- which(apply(dat, 2, function(x) all(is.na(x)))) - 1
  m <- ncol(dat) - 1
  missing_tail <- seq(1 + m - length(all_missing), m)
  orig_order <- seq_len(m + 1) - 1
  new_order <- c(orig_order[!orig_order %in% all_missing], all_missing)
  revert_order <- order(new_order) - 1

  missing_vec <- apply(dat, 1, function(x) which(is.na(x)) - 1)
  missing_vec_unique <- unique(missing_vec)
  missing_pattern <- sapply(missing_vec, paste, collapse = '_')
  missing_pattern_unique <- unique(missing_pattern)
  pattern_inds <- lapply(missing_pattern_unique, function(x) which(missing_pattern == x) - 1)
  present_vec <- apply(dat, 1, function(x) which(!is.na(x)) - 1)
  present_vec_unique <- unique(present_vec)
  n_pattern <- length(missing_pattern_unique)

  list(
    # Number of different patterns
    npatt = n_pattern,
    # Distinct patterns of missing-ness
    mlist = missing_vec_unique,
    # Distinct patterns of present-ness
    plist = present_vec_unique,
    # Indices matching patterns to rows
    indlist = pattern_inds,
    # Columns where all values are missing
    all_missing = all_missing,
    # Stick all of the "all_missing" data at the end
    missing_tail = missing_tail,
    # New order, with all-missing data at the end
    new_order = new_order,
    # Indices for reverting back to the original order
    revert_order = revert_order
  )

}
