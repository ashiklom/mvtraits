rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

global_mean <- c(-5, 0, 5)
nparam <- length(global_mean)
global_cor <- matrix(nrow = 3, ncol = 3,
                     data = c(1, 0.75, 0.25,
                              0.75, 1, 0,
                              0.25, 0, 1))
global_sigmas <- diag(c(1, 2, 3))
global_cov <- global_sigmas %*% global_cor %*% global_sigmas

ngroup <- 10
group_mu <- random_mvnorm(ngroup, global_mean, global_cov)
group_sigma <- diag(nparam)

group_cor_estimate <- cor(group_mu)
print(group_cor_estimate)

ndat <- 5000
dat <- matrix(0, ndat, nparam)
groups <- sample.int(ngroup, ndat, TRUE)
for (i in seq_len(ndat)) {
    dat[i,] <- random_mvnorm(1, group_mu[groups[i],], group_sigma)
}

niter <- 500
nchains <- 2
parallel <- TRUE

message('Running simple hierarchical...')
samps_hier <- fit_mvnorm_hier(dat, groups, niter = niter, nchains = nchains, parallel = parallel,
                              autofit = TRUE)
message('Done!')

samps_hier_full <- add_correlations(samps_hier, hier = TRUE, ngroups = ngroup)
samps_hier_mcmc <- results2mcmclist(samps_hier_full, 'hier')
samps_hier_burned <- window(samps_hier_mcmc, start = floor(niter / 2))
hier_sum <- summary_df(samps_hier_burned, group = TRUE)

if (interactive()) {
    library(ggplot2)
    ggplot(dplyr::filter(hier_sum, variable == 'Corr')) + 
        aes(x = group, y = Mean, ymin = `2.5%`, ymax = `97.5%`, color = group) + 
        geom_pointrange() +
        facet_wrap(~index)
}
