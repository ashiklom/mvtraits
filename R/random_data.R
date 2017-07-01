#' @export
random_data <- function(mu = colMeans(iris[,-5]), Sigma = cov(iris[,-5]),
                        N = 5000, ngroup = 7, frac_miss = 0.75) {
    Cor <- cov2cor(Sigma)
    dat_all <- mvtnorm::rmvnorm(N, mu, Sigma)
    groups <- sample.int(ngroup, N, replace = TRUE)
    dat <- dat_all
    nmiss <- round(length(dat) * frac_miss)
    miss <- sample.int(length(dat), size = nmiss)
    dat[miss] <- NA
    out <- list(dat = dat, dat_all = dat_all, groups = groups,
                mu = mu, Sigma = Sigma, Cor = Cor)
    return(out)
}
