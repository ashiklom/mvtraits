library(Amelia)
library(tidyverse)
library(magrittr)

trydat <- readRDS('~/Projects/new-phytologist-traits/np_trait_analysis/traits_analysis.rds')

select_rxp <- function(area_mass) {
    area_rxp <- 'leaf_lifespan|LMA|area'
    mass_rxp <- 'leaf_lifespan|SLA|mass'
    use_rxp <- switch(area_mass, area = area_rxp, mass = mass_rxp)
    return(use_rxp)
}

get_df <- function(dat, area_mass) {
    rxp <- select_rxp(area_mass)
    dat %>% 
        dplyr::select(pft, dplyr::matches(rxp)) %>% 
        #dplyr::select(-dplyr::matches('Jmax|Vcmax|Rd')) %>% 
        dplyr::filter_at(dplyr::vars(-pft), dplyr::any_vars(!is.na(.)))
}

get_nested <- function(dat, area_mass) {
    get_df(dat, area_mass) %>% tidyr::nest(-pft)
}

get_datamatrix <- function(dat, area_mass) {
    data_df <- get_df(dat, area_mass)
    data_mat <- data_df %>% dplyr::select(-pft) %>% as.matrix() %>% log10()
    data_groups <- data_df %>% dplyr::pull(pft) %>% as.integer()
    return(list(dat = data_mat, groups = data_groups))
}

area_nest <- get_nested(trydat, 'area')

area_imp <- area_nest %>% 
    mutate(imputed = map(data, amelia))
dir.create('data')
saveRDS(area_imp, 'data/amelia_pft.rds')
