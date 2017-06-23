#library(mvtraits)
devtools::load_all('.')
library(tidyverse)

options(digits = 3)

trydat <- readRDS('~/Projects/new-phytologist-traits/np_trait_analysis/traits_analysis.rds')

get_datamatrix <- function(dat, area_mass) {
    area_rxp <- 'leaf_lifespan|LMA|area'
    mass_rxp <- 'leaf_lifespan|SLA|mass'
    use_rxp <- switch(area_mass, area = area_rxp, mass = mass_rxp)
    data_df <- dat %>% 
        dplyr::select(pft, dplyr::matches(use_rxp)) %>% 
        #dplyr::select(-dplyr::matches('Jmax|Vcmax|Rd')) %>% 
        dplyr::filter_at(dplyr::vars(-pft), dplyr::any_vars(!is.na(.)))
    data_mat <- data_df %>% dplyr::select(-pft) %>% as.matrix() %>% log10()
    data_groups <- data_df %>% dplyr::pull(pft) %>% as.integer()
    return(list(dat = data_mat, groups = data_groups))
}

niter <- 500
nchain <- 1
parallel <- FALSE

area_data <- get_datamatrix(trydat, 'area')
dat <- area_data[['dat']]
groups <- area_data[['groups']]
ngroup <- length(unique(groups))
dat_sub <- dat

nparam <- ncol(dat_sub)
strong_wish <- Wishart_prior_param(1, 1/10000, nparam)
priors <- list(v_global = strong_wish$v0, S_global = strong_wish$S0,
               Sigma_global = diag(1/1000, nparam))

Sigma0_inv <- solve(priors[['Sigma_global']])

area_fit_multivariate <- fit_mvnorm(dat_sub, priors = priors, 
                                    niter = niter, nchains = nchain, parallel = parallel)

#area_fit_mcmc <- area_fit_multivariate %>% 
    #add_correlations() %>% 
    #results2mcmclist(chain2matrix_multi) %>% 
    #window(start = floor(niter / 2))
#plot(area_fit_mcmc, ask = TRUE)

#priors[['v_group']] <- rep(strong_wish$v0, ngroup)
#priors[['S_group']] <- matrep(strong_wish$S0, ngroup)

#t1 <- proc.time()
#area_fit_hierarchical <- fit_mvnorm_hier(dat_sub, groups = groups[groups %in% 1:4], 
                                         #niter = niter, nchains = nchain, parallel = parallel)
#t2 <- proc.time()
#print(t2 - t1)

#area_fit_mcmc <- area_fit_hierarchical %>% 
    #add_correlations() %>% 
    #results2mcmclist(chain2matrix_hier)
#plot(area_fit_mcmc, ask = TRUE)

#area_mv_mean <- apply(area_fit_multivariate$mu, )

#mass_fit_multivariate <- fit_mvnorm(mass_data[['dat']], niter = 10000, nchains = 1, parallel = FALSE)
#mass_fit_mcmc <- mass_fit_multivariate %>%
    #add_correlations() %>% 
    #results2mcmclist(chain2matrix_multi) %>% 
    #window(start = floor(niter / 2))

#area_fit_hier <- fit_mvnorm_hier(area_data[['dat']], area_data[['groups']])
#mass_fit_hier <- fit_mvnorm_hier(mass_data[['dat']], mass_data[['groups']])
