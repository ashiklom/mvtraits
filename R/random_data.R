#' Random simple multivariate data
#'
#' Useful for testing.
#' @param mu Vector of means
#' @param sigma Covariance matrix
#' @param n Number of samples
#' @param all_missing Column containing all missing data. Ignored if `NA`.
#' @inheritParams random_missing
#' @export
random_data_multi <- function(mu = c(-2, 0, 2),
                              sigma = matrix(c(1, 0.9, 0.6,
                                               0.9, 1, 0.3,
                                               0.6, 0.3, 1), 3, 3),
                              n = 5000,
                              frac_miss = 0.75,
                              all_missing = 3) {
    cor <- cov2cor(sigma)
    dat_all <- random_mvnorm(n, mu, sigma)
    dat <- random_missing(dat_all, frac_miss)
    if (is.numeric(all_missing)) {
        dat[, all_missing] <- NA
    }
    list(dat = dat, dat_all = dat_all,
         mu = mu, sigma = sigma, cor = cor)
}

# Set up defaults for data generation
# Math for the bivariate case is easier
def_mu_global <- c(x = 0, y = 0)
def_sigma_global <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
dimnames(def_sigma_global) <- rep(list(names(def_mu_global)), 2)
slope <- def_sigma_global[2, 1] * def_sigma_global[2, 2] / def_sigma_global[1, 1]
mu_group_x <- c(-0.75, 0, 0.75)
def_mu_group <- as.list(data.frame(rbind(mu_group_x, slope * mu_group_x)))
mu_group_mat <- Reduce(rbind, def_mu_group)

sd_group <- c(0.2, 0.2)

cor_group_same <- lapply(
    list(a = def_sigma_global, b = def_sigma_global, c = def_sigma_global),
    cov2cor
)

cor_group_diff <- list(
    a = matrix(c(1, 0.75, 0.75, 1), nrow = 2),
    b = matrix(c(1, 0, 0, 1), nrow = 2),
    c = matrix(c(1, -0.75, -0.75, 1), nrow = 2)
)

cor2cov <- function(cor, sdvec) {
    sddiag <- diag(sdvec)
    sddiag %*% cor %*% sddiag
}

cov_group_same <- lapply(cor_group_same, cor2cov, sdvec = sd_group)
def_sigma_group <- lapply(cor_group_diff, cor2cov, sdvec = sd_group)

#' Random hierarchical multivariate data
#'
#' Useful for testing.
#' @param mu_global Global mean vector
#' @param sigma_global Global covariance matrix
#' @param mu_group Group mean vectors (list)
#' @param sigma_group Group covariance matrices (list)
#' @param all_missing_group Group containing a completely missing column
#' @param all_missing_col Column that is missing completely
#' @inheritParams random_missing
#' @export
random_data_hier <- function(mu_global = def_mu_global,
                             sigma_global = def_sigma_global,
                             mu_group = def_mu_group,
                             sigma_group = def_sigma_group,
                             n = 5000,
                             frac_miss = 0.75,
                             all_missing_group = 3,
                             all_missing_col = 2
) {
    ngroups <- length(mu_group)
    groups <- rep(seq_len(ngroups), each = floor(n / ngroups), length.out = n)
    samps_list <- Map(
        random_mvnorm,
        mu = mu_group,
        Sigma = sigma_group,
        n = ceiling(n / ngroups)
    )
    stopifnot(
      all_missing_group %in% groups,
      length(groups) == n
    )
    dat_all <- do.call(rbind, samps_list)[seq_len(n), ]
    dat <- random_missing(dat_all, frac_miss)
    dat[groups == all_missing_group, all_missing_col] <- NA
    list(dat = dat, dat_all = dat_all, groups = groups,
         mu_global = mu_global, sigma_global = sigma_global,
         cor_global = cov2cor(sigma_global),
         mu_group = mu_group, sigma_group = sigma_group,
         cor_group = lapply(sigma_group, cov2cor)
    )
}

#' Intersperse a random number of NAs through data
#'
#' @param dat Data matrix to which to add `NA`
#' @param frac_miss Fraction of randomly missing data
#' @export
random_missing <- function(dat, frac_miss = 0.75) {
  stopifnot(is.matrix(dat))
  nmiss <- floor(length(dat) * frac_miss)
  miss <- sample.int(length(dat), size = nmiss)
  dat[miss] <- NA
  dat
}
