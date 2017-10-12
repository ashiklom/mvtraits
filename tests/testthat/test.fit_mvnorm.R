rm(list = ls())

library(testthat)
library(mvtraits)

rand <- random_data()

niter <- 10000
nchains <- 2
parallel <- FALSE
keep_samples <- Inf
max_attempts <- 20
autofit <- TRUE

message('Running simple multivariate...')
# debugonce(fit_mvnorm)
rand$dat[,4] <- NA
result <- fit_mvnorm(rand$dat, niter = niter, nchains = nchains, parallel = parallel,
  autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts)
# R version: 25 seconds
# R main, C functions: 11 seconds
# All C:

message('Done!')
