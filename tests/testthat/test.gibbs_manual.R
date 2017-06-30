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
Cor <- cov2cor(Sig)

N <- 5000
ngroup <- 7
frac_miss <- 0.9
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)
groups <- sample.int(ngroup, N, replace = TRUE)

# Randomly remove a fraction of the data
dat <- dat_all

nmiss <- round(length(dat) * frac_miss)
miss <- sample.int(length(dat), size = nmiss)

dat[miss] <- NA

niter <- 2000
nchain <- 1
parallel <- FALSE

message('Running simple multivariate...')
samps_mv <- fit_mvnorm(dat, niter = niter, nchains = nchain, parallel = parallel)
message('Done!')

samps_mv_full <- add_correlations(samps_mv)
samps_mv_mcmc <- results2mcmclist(samps_mv_full, chain2matrix_multi)
samps_mv_burned <- window(samps_mv_mcmc, start = floor(niter / 2))
mv_sum <- summary_df(samps_mv_burned, group = NULL)

message('Running hierarchical...')
samps_hier <- fit_mvnorm_hier(dat, groups, niter = niter, nchains = nchain, parallel = parallel)
message('Done!')

samps_hier_full <- add_correlations(samps_hier)
samps_hier_mcmc <- results2mcmclist(samps_hier_full, chain2matrix_hier)
samps_hier_burned <- window(samps_hier_mcmc, start = floor(niter / 2))
hier_sum <- summary_df(samps_hier_burned, group = TRUE)

library(ggplot2)
ggplot(dplyr::filter(hier_sum, variable == 'Corr')) + 
    aes(x = group, y = Mean, ymin = `2.5%`, ymax = `97.5%`, color = group) + 
    geom_pointrange() +
    facet_wrap(~index, scales = 'free_y')
