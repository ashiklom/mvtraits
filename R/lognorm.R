#' Lognormal distribution parameter conversion functions
#'
#' Calculate lognormal distribution mean and variance from its log-scale 
#' parameters.
#'
#' @param mu Log-scale mean
#' @param sigma2 Log-scale variance
#' @param base Logarithm base. Default is `e` (`exp(1)`)
#' @export
log_mean <- function(mu, sigma2, base = exp(1)) {
  base ^ (mu + sigma2 / 2)
}

#' @rdname log_mean
#' @export
log_var <- function(mu, sigma2, base = exp(1)) {
  (base ^ sigma2 - 1) * base ^ (2 * mu + sigma2)
}

#' @rdname log_mean
#' @export
log10_mean <- function(mu, sigma2) {
  log_mean(mu, sigma2, base = 10)
}

#' @rdname log_mean
#' @export
log10_var <- function(mu, sigma2) {
  log_var(mu, sigma2, base = 10)
}
