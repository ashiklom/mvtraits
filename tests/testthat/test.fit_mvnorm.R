rm(list = ls())

library(testthat)
library(mvtraits)

rand <- random_data()

niter <- 10000
nchains <- 2
parallel <- FALSE
keep_samples <- Inf
max_attempts <- 20

message('Running simple multivariate...')
# debugonce(fit_mvnorm)
rand$dat[,4] <- NA
t1 <- proc.time()
samps_mv <- fit_mvnorm(rand$dat, niter = niter, nchains = nchains, parallel = parallel,
                       autofit = FALSE, keep_samples = keep_samples, max_attempts = max_attempts)
t2 <- proc.time()
print(t2 - t1)
# R version: 25 seconds
# R main, C functions: 11 seconds
# All C:

samps_mcmc <- results2mcmclist(samps_mv, 'multi')
plot(samps_mcmc[,"Sigma..par04..par04"])

message('Done!')


samps_mv_full <- add_correlations(samps_mv)
samps_mv_mcmc <- results2mcmclist(samps_mv_full, 'multi')
samps_mv_burned <- window(samps_mv_mcmc, start = floor(niter / 2))
mv_sum <- summary_df(samps_mv_burned, group = NULL)

if (exists('doplot')) {
    plot(samps_mv_mcmc, ask = TRUE)
}
