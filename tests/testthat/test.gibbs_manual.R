rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

# Simulate some data
mu <- colMeans(iris[,-5])
Sig <- cov(iris[,-5])

N <- 5000
ngroup <- 7
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)
groups <- sample.int(ngroup, N, replace = TRUE)

# Randomly remove a fraction of the data
dat <- dat_all

nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)

dat[miss] <- NA

message('Running simple multivariate...')
niter <- 5000
samps_mv <- fit_mvnorm(dat, niter = niter)
message('Done!')

samps_mv_mcmc <- results2mcmclist(samps_mv, chain2matrix_multi)
samps_mv_burned <- window(samps_mv_mcmc, start = floor(niter / 2))
summary(samps_mv_burned)

message('Running hierarchical...')
samps_hier <- fit_mvnorm_hier(dat, groups)
message('Done!')

samps_hier_mcmc <- results2mcmclist(samps_hier, chain2matrix_hier)
samps_hier_burned <- window(samps_hier_mcmc, start = floor(niter / 2))
summary(samps_hier_burned)
