rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

rand <- random_data()

niter <- 500
nchains <- 2
parallel <- TRUE

message('Running simple hierarchical...')
samps_hier <- fit_mvnorm_hier(rand$dat, rand$groups, niter = niter, nchains = nchains, parallel = parallel,
                              repeat_until_converged = TRUE)
message('Done!')

samps_hier_full <- add_correlations(samps_hier)
samps_hier_mcmc <- results2mcmclist(samps_hier_full, chain2matrix_hier)
samps_hier_burned <- window(samps_hier_mcmc, start = floor(niter / 2))
hier_sum <- summary_df(samps_hier_burned, group = TRUE)

if (interactive()) {
    library(ggplot2)
    ggplot(dplyr::filter(hier_sum, variable == 'mu')) + 
        aes(x = group, y = Mean, ymin = `2.5%`, ymax = `97.5%`, color = group) + 
        geom_pointrange() +
        facet_wrap(~index, scales = 'free_y')
}
