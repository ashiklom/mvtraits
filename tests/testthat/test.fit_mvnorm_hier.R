rm(list = ls())

library(testthat)
library(mvtraits)
plt <- FALSE

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

niter <- 5000
nchains <- 2
parallel <- FALSE
autofit <- TRUE
keep_samples <- 2000

message('Running simple hierarchical...')
result <- fit_mvnorm_hier(dat, groups, niter = niter, nchains = nchains, parallel = parallel,
                              autofit = autofit, keep_samples = keep_samples)
message('Done!')


#if (exists('doplot')) {
    library(ggplot2)
    ggplot(dplyr::filter(hier_sum, variable == 'Corr')) +
        aes(x = group, y = Mean, ymin = `2.5%`, ymax = `97.5%`, color = group) +
        geom_pointrange() +
        facet_wrap(~index)
#}
