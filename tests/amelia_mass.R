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

mass_nest <- get_nested(trydat, 'mass')

mass_imp <- mass_nest %>% 
    mutate(imputed = map(data, amelia))
dir.create('data')
saveRDS(mass_imp, 'data/amelia_pft_mass.rds')

arraycorr <- function(arr) {
    out <- arr
    for (i in seq_len(dim(arr)[3])) {
        out[,,i] <- cov2cor(arr[,,i])
    }
    return(out)
}

mass_imp2 <- mass_imp %>% 
    filter(sapply(imputed, class) == 'amelia') %>% 
    mutate(mu_all = map(imputed, 'mu'),
           mu_all = map2(mu_all, data, ~"rownames<-"(.x, colnames(.y))),
           mu_means = map(mu_all, rowMeans),
           Sigma_all = map(imputed, 'covMatrices'),
           Sigma_all = map2(Sigma_all, data, ~"dimnames<-"(.x, list(colnames(.y), colnames(.y), NULL))),
           Sigma_means = map(Sigma_all, apply, 1:2, mean),
           Corr_all = map(Sigma_all, arraycorr),
           Corr_all = map2(Corr_all, data, ~"dimnames<-"(.x, list(colnames(.y), colnames(.y), NULL))),
           Corr_means = map(Corr_all, apply, 1:2, mean))

vec2df <- function(vec, ...) tibble(index = names(vec), value = vec, ...)
mat2df <- function(mat, diag = TRUE, ...) {
    pars <- rownames(mat)
    lt <- which(lower.tri(mat, diag = diag), arr.ind = TRUE)
    matvec <- mat[lt]
    names(matvec) <- paste(pars[lt[,1]], pars[lt[,2]], sep = '..')
    vec2df(matvec, ...)
}

mass_proc <- mass_imp2 %>% 
    mutate(mu_df = map(mu_means, vec2df, variable = 'mu'),
           sigma_df = map(Sigma_means, mat2df, variable = 'sigma'),
           corr_df = map(Corr_means, mat2df, variable = 'corr', diag = FALSE)) %>% 
    mutate(alldat = pmap(list(mu_df, sigma_df, corr_df), bind_rows)) %>% 
    select(pft, alldat) %>% 
    unnest()

mass_means_wide_all <- mass_proc %>% 
    filter(variable == 'mu') %>% 
    select(-variable) %>% 
    spread(index, value)

mass_means_wide <- select(mass_means_wide_all, -pft)

mass_glob_mean <- vec2df(colMeans(mass_means_wide), pft = 'global', variable = 'mu')
mass_glob_cov <- mat2df(cov(mass_means_wide), pft = 'global', variable = 'sigma')
mass_glob_cor <- mat2df(cor(mass_means_wide), diag = FALSE, pft = 'global', variable = 'corr')
mass_glob <- bind_rows(mass_glob_mean, mass_glob_cov, mass_glob_cor)
mass_all <- full_join(mass_glob, mass_proc) %>% 
    mutate(pft = factor(pft, levels = c('global', levels(mass_proc[['pft']])))) %>% 
    separate(index, c('xvar', 'yvar'), sep = '\\.\\.', remove = FALSE)

npft <- n_distinct(mass_all$pft)
pfts <- levels(mass_all$pft)
pfts <- pfts[pfts %in% mass_all$pft]
pft_colors <- c('black', RColorBrewer::brewer.pal(npft - 1, 'Set1'))
names(pft_colors) <- pfts

############################################################
# Bar chart
dir.create('figures')
pdf('figures/mass_corrmat.pdf', width = 12, height = 10)
ggplot(mass_all %>% filter(variable == 'corr')) + 
    aes(x = pft, y = value, fill = pft) + 
    geom_col() + 
    facet_wrap(~xvar + yvar, drop = FALSE) + 
    scale_fill_manual(values = pft_colors) + 
    theme(axis.text.x = element_blank()) + 
    xlab('Plant functional type') + 
    ylab('Mean correlation coefficient')
dev.off()

############################################################
# Ellipse figure
mass_mean <- mass_all %>% 
    filter(variable == 'mu') %>% 
    select(xvar, value, pft) %>% 
    spread(xvar, value)

mass_both <- mass_imp2 %>% 
    select(pft, mu_means, Sigma_means) %>% 
    bind_rows(tibble(pft = 'global', 
                     mu_means = list(colMeans(mass_means_wide)),
                     Sigma_means = list(cov(mass_means_wide))))

ellipsePlot <- function(xvar, yvar, radius = qnorm(0.95)) {
    vrs <- c(xvar, yvar)
    plt <- mass_both %>% 
        mutate(mus = map(mu_means, ~.[vrs]),
               mu_df = map(mus, ~tibble(xvar = .x[xvar], yvar = .x[yvar])),
               sigmas = map(Sigma_means, ~.[vrs, vrs]),
               ellipses = pmap(list(center = mus, shape = sigmas), car::ellipse, radius = radius, 
                               draw = FALSE),
               ellipses = map(ellipses, as_data_frame))
    gplot <- ggplot() + 
        aes(color = pft) + 
        geom_path(data = unnest(plt, ellipses), aes(x = x, y = y)) + 
        geom_point(data = unnest(plt, mu_df), aes(x = xvar, y = yvar)) +
        xlab(xvar) + 
        ylab(yvar) + 
        scale_color_manual(values = pft_colors)
    return(gplot)
}

trait_combinations <- distinct(mass_all %>% filter(variable == 'corr'), xvar, yvar)

plots <- trait_combinations %$% 
    map2(xvar, yvar, ellipsePlot, radius = qnorm(0.9))
files <- trait_combinations %$%
    file.path('figures', paste(xvar, yvar, 'pdf', sep = '.'))
walk2(files, plots, ggsave, device = pdf, width = 10, height = 8)
