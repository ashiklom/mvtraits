#' Bootstrap missing values in data from fitted multivariate model
#'
#' @param fit Fitting result, as returned by [fit_mvnorm()] or [fit_mvnorm_hier()].
#' @param dat Data with partially missing values
#' @param n Number of bootstrap samples (default = 500)
#' @return If `missing_only`, matrix of bootstrapped missing values (`nmissing x
#'   n`). Otherwise, an array of bootstraps of the full imputed dataset, with
#'   dimensions `c(nrow(dat), ncol(dat), n)`.
#' @author Alexey Shiklomanov
#' @export
bootstrap_missing <- function(fit, dat, n = 500) {
  samp <- fit$samples
  stopifnot(
    inherits(samp, "mcmc.list"),
    any(is.na(dat))
  )
  if (!is.matrix(dat)) dat <- as.matrix(dat)
  samp_mat <- as.matrix(samp)
  # Grab n samples from the posterior
  fit_samp <- samp_mat[sample(nrow(samp_mat), n), ]
  fit_mu <- fit_samp[, grep("^mu\\.\\.", colnames(fit_samp))]
  fit_sigma_wide <- fit_samp[, grep("^Sigma\\.\\.", colnames(fit_samp))]
  setup <- setup_missing(dat)
  out <- array(numeric(), c(nrow(dat), ncol(dat), n))
  for (i in seq_len(n)) {
    mu <- fit_mu[i, ]
    sigma <- lowerdiag2mat(fit_sigma_wide[i, ], col_names = FALSE, corr = FALSE,
                           hier = hier)
    fill <- mvnorm_fill_missing(dat, mu, sigma, setup = setup)
    out[, , i] <- fill
  }
  out
}

#' Bootstrap missing values in data from fitted multivariate hierarchical model
#'
#' @param groups a vector of groups. If `NULL`
#'   (default), assume a simple multivariate fit.
#' @author Alexey Shiklomanov
#' @inheritParams bootstrap_missing
#' @inherit bootstrap_missing return
#' @export
bootstrap_missing_hier <- function(fit, dat, groups, n = 500) {
  samp <- fit$samples
  stopifnot(
    inherits(samp, "mcmc.list"),
    any(is.na(dat)),
    length(groups) == nrow(dat)
  )
  if (!is.matrix(dat)) dat <- as.matrix(dat)
  samp_mat <- as.matrix(samp)
  ugroup <- unique(groups)
  # Grab n samples from the posterior
  fit_samp <- samp_mat[sample(nrow(samp_mat), n), ]
  fit_mu_l <- lapply(ugroup, get_mu_group, smat = fit_samp)
  names(fit_mu_l) <- ugroup
  fit_sigma_wide_l <- lapply(ugroup, get_sigma_group, smat = fit_samp)
  names(fit_sigma_wide_l) <- ugroup
  setup_l <- lapply(
    ugroup,
    function(x) setup_missing(dat[groups == x, ])
  )
  out <- array(numeric(), c(nrow(dat), ncol(dat), n))
  for (i in seq_len(n)) {
    for (g in ugroup) {
      g_mu <- fit_mu_l[[g]][i,]
      g_sigma <- lowerdiag2mat(
        fit_sigma_wide_l[[g]][i, ],
        col_names = FALSE,
        corr = FALSE,
        hier = hier
      )
      g_rows <- which(groups == g)
      g_dat <- dat[g_rows, ]
      fill <- mvnorm_fill_missing(g_dat, g_mu, g_sigma, setup = setup_l[[g]])
      out[g_rows, , i] <- fill
    }
  }
  out
}

get_mu_group <- function(group, smat) {
  rx <- paste0("^mu\\.\\.", group, "\\.\\.")
  i <- grep(rx, colnames(smat))
  smat[, i]
}

get_sigma_group <- function(group, smat) {
  rx <- paste0("^Sigma\\.\\.", group, "\\.\\.")
  i <- grep(rx, colnames(smat))
  smat[, i]
}
