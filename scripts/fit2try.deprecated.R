library(mvtraits)
#devtools::load_all('.')
library(tidyverse)

options(digits = 3)

trydat <- readRDS('extdata/traits_analysis.rds')

results_dir <- 'results'
dir.create(results_dir, showWarnings = FALSE)

progress_dir <- 'progress'
dir.create(progress_dir, showWarnings = FALSE)

get_datamatrix <- function(dat, area_mass) {
    area_rxp <- 'leaf_lifespan|SLA|area'
    mass_rxp <- 'leaf_lifespan|SLA|mass'
    use_rxp <- switch(area_mass, area = area_rxp, mass = mass_rxp)
    data_df <- dat %>% 
        dplyr::select(pft, dplyr::matches(use_rxp)) %>% 
        dplyr::filter_at(dplyr::vars(-pft), dplyr::any_vars(!is.na(.)))
    data_mat <- data_df %>% dplyr::select(-pft) %>% as.matrix() %>% log10()
    data_groups <- data_df %>% dplyr::pull(pft) %>% as.integer()
    return(list(dat = data_mat, groups = data_groups))
}

niter <- 1000
nchain <- 5
parallel <- TRUE

mass_data <- get_datamatrix(trydat, 'area')
dat <- mass_data[['dat']]
groups <- mass_data[['groups']]
ngroup <- length(unique(groups))
dat_sub <- dat

#nparam <- ncol(dat_sub)
#strong_wish <- Wishart_prior_param(1, 1/10000, nparam)
#priors <- list(v_global = strong_wish$v0, S_global = strong_wish$S0,
               #Sigma_global = diag(1/1000, nparam))
priors <- list()

#message('Starting multivariate simulation')
#area_fit_multivariate <- fit_mvnorm(dat_sub, priors = priors, 
                                    #niter = niter, nchains = nchain, parallel = parallel,
                                    #autofit = TRUE, 
                                    #save_progress = file.path(progress_dir, 'area_multi_global'))
#message('Done!')

#save_path <- file.path(results_dir, 'area_multi_global.rds')
#saveRDS(area_fit_multivariate, save_path)

#priors[['v_group']] <- rep(strong_wish$v0, ngroup)
#priors[['S_group']] <- matrep(strong_wish$S0, ngroup)

message('Starting hierarchical sampling')
mass_fit_hierarchical <- fit_mvnorm_hier(dat, groups = groups, niter = niter,
                                         nchains = nchain, parallel = parallel,
                                         autofit = TRUE, max_attempts = 200, 
                                         save_progress = file.path(progress_dir, 'area_hier'))
message('Done!')
saveRDS(mass_fit_hierarchical, file.path(results_dir, 'area_hier.rds'))

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
