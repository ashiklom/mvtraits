rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

rand <- random_data()

niter <- 2000
nchains <- 1
parallel <- FALSE

message('Running simple multivariate...')
samps_mv <- fit_mvnorm(rand$dat, niter = niter, nchains = nchains, parallel = parallel)
message('Done!')

samps_mv_full <- add_correlations(samps_mv)
samps_mv_mcmc <- results2mcmclist(samps_mv_full, chain2matrix_multi)
samps_mv_burned <- window(samps_mv_mcmc, start = floor(niter / 2))
mv_sum <- summary_df(samps_mv_burned, group = NULL)

if (interactive()) {
    plot(samps_mv_mcmc, ask = TRUE)
}
