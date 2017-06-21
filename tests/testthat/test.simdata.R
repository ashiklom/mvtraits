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

N <- 100
ngroup <- 3
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)
groups <- sample.int(ngroup, N, replace = TRUE)

# Randomly remove half of the data
dat <- dat_all

nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)

dat[miss] <- NA

custom_inputs <- list()

#fit_uni <- runModel('uni', dat[groups == 1,], iter = 100, max.attempts = 1)
#fit_multi <- runModel('multi', dat[groups == 1,], iter = 100, max.attempts = 1)
#fit_hier <- runModel('hier', dat, groups = groups, iter = 100, max.attempts = 1)

file.remove(list.files('testmodel_*.rds'))
