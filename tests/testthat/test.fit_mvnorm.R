rm(list = ls())

library(testthat)
library(mvtraits)
plt <- FALSE

rand <- random_data()

niter <- 100
nchains <- 2
parallel <- FALSE

message('Running simple multivariate...')
samps_mv <- fit_mvnorm(rand$dat, niter = niter, nchains = nchains, parallel = parallel, 
                       autofit = TRUE)
message('Done!')

samps_mv_full <- add_correlations(samps_mv)
samps_mv_mcmc <- results2mcmclist(samps_mv_full, 'multi')
samps_mv_burned <- window(samps_mv_mcmc, start = floor(niter / 2))
mv_sum <- summary_df(samps_mv_burned, group = NULL)

if (exists('doplot')) {
    plot(samps_mv_mcmc, ask = TRUE)
}
