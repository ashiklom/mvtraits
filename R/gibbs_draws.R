#' Gibbs draw of mean vector
#'
#' @param data Matrix of data values
#' @param Sigma Current draw of covariance matrix
#' @param mu0 Prior mean values
#' @param Sigma0 Prior mean covariance
#' @export
draw_mu <- function(data, Sigma, mu0, Sigma0) {
  Sigma_inv <- solve(Sigma)
  Sigma_0_inv <- solve(Sigma_0)
  xbar <- colMeans(data)
  nx <- nrow(data)
  c_draw_mu(xbar, nx, Sigma_inv, mu0, Sigma0)
}

#' Gibbs draw of covariance matrix
#'
#' @param data Matrix of data values
#' @param mu Current draw of mean vector
#' @param v0 Prior Wishart degrees of freedom
#' @param S0 Prior Wishart rate matrix
#' @export
draw_Sigma <- function(data, mu, v0, S0) {
  c_draw_Sigma(data, mu, v0, S0)
}
