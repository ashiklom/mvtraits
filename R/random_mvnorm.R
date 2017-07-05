#' @export
random_mvnorm <- function(n, mu, Sigma) {
    if (!is.matrix(mu)) {
        mu <- t(mu)
    }
    c_random_mvnorm(n, mu, Sigma)
}
