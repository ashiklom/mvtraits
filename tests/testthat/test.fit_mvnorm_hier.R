rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

rand <- random_data(ngroup = 3)
attach(rand)
ngroups <- length(unique(groups))

niter <- 500
nchains <- 2
parallel <- FALSE

message('Running simple hierarchical...')
samps_hier <- fit_mvnorm_hier(dat, groups, niter = niter, nchains = nchains, parallel = parallel,
                              autofit = TRUE)
message('Done!')

samps_hier_full <- add_correlations(samps_hier, hier = TRUE, ngroups = ngroups)
samps_hier_mcmc <- results2mcmclist(samps_hier_full, 'hier')
samps_hier_burned <- window(samps_hier_mcmc, start = floor(niter / 2))
hier_sum <- summary_df(samps_hier_burned, group = TRUE)

if (interactive()) {
    library(ggplot2)
    ggplot(dplyr::filter(hier_sum, variable == 'mu')) + 
        aes(x = group, y = Mean, ymin = `2.5%`, ymax = `97.5%`, color = group) + 
        geom_pointrange() +
        facet_wrap(~index, scales = 'free_y')
}
