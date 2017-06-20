rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

# Simulate some data
mu <- c(10, 5, 0, -5, 10)
Sig <- clusterGeneration::genPositiveDefMat(length(mu))$Sigma

N <- 1000
ngroup <- 7
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)
groups <- sample.int(ngroup, N, replace = TRUE)

# Randomly remove half of the data
dat <- dat_all

nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)

dat[miss] <- NA

samps_mv <- fit_mvnorm(dat)
samps_hier <- fit_mvnorm_hier(dat, groups)
