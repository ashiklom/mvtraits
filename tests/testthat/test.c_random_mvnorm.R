library(testthat)
devtools::load_all('.')

ss <- function(x, y) {
    sum((y - x)^2)
}

mu <- colMeans(iris[,-5])
Sigma <- cov(iris[,-5])

n <- 5000
m <- length(mu)

samp <- mvtraits:::c_random_mvnorm(n, t(mu), Sigma)
samp_mu <- colMeans(samp)
samp_cov <- cov(samp)

test_that('c_random_mvnorm (vector mu) produces reasonable values', {
             expect_lt(ss(samp_mu, mu), 0.01)
             expect_lt(ss(samp_cov, Sigma), 0.05)
})

mu_mat <- matrix(mu, n, m, byrow = TRUE)
samp2 <- mvtraits:::c_random_mvnorm(n, mu_mat, Sigma)
samp2_mu <- colMeans(samp2)
samp2_cov <- cov(samp2)

test_that('c_random_mvnorm (matrix mu) produces reasonable values', {
             expect_lt(ss(samp2_mu, mu), 0.01)
             expect_lt(ss(samp2_cov, Sigma), 0.05)
})
