library(Amelia)
library(tidyverse)
library(magrittr)

area_imp <- readRDS('data/amelia_pft.rds')

arraycorr <- function(arr) {
    out <- arr
    for (i in seq_len(dim(arr)[3])) {
        out[,,i] <- cov2cor(arr[,,i])
    }
    return(out)
}

trait_order <- c('leaf_lifespan', 'LMA', 'Narea', 'Parea', 'Rdarea', 'Vcmax_area', 'Jmax_area')

area_imp2 <- area_imp %>% 
    filter(sapply(imputed, class) == 'amelia') %>% 
    mutate(mu_all = map(imputed, 'mu'),
           mu_all = map2(mu_all, data, ~"rownames<-"(.x, colnames(.y))),
           mu_means = map(mu_all, rowMeans),
           Sigma_all = map(imputed, 'covMatrices'),
           Sigma_all = map2(Sigma_all, data, ~"dimnames<-"(.x, list(colnames(.y), colnames(.y), NULL))),
           Sigma_all = map(Sigma_all, ~.x[trait_order,trait_order,]),
           Sigma_means = map(Sigma_all, apply, 1:2, mean),
           Corr_all = map(Sigma_all, arraycorr),
           Corr_means = map(Corr_all, apply, 1:2, mean))

vec2df <- function(vec, ...) tibble(index = names(vec), value = vec, ...)
mat2df <- function(mat, diag = TRUE, ...) {
    pars <- rownames(mat)
    lt <- which(lower.tri(mat, diag = diag), arr.ind = TRUE)
    matvec <- mat[lt]
    names(matvec) <- paste(pars[lt[,1]], pars[lt[,2]], sep = '..')
    vec2df(matvec, ...)
}

area_proc <- area_imp2 %>% 
    mutate(mu_df = map(mu_means, vec2df, variable = 'mu'),
           sigma_df = map(Sigma_means, mat2df, variable = 'sigma'),
           corr_df = map(Corr_means, mat2df, variable = 'corr', diag = FALSE)) %>% 
    mutate(alldat = pmap(list(mu_df, sigma_df, corr_df), bind_rows)) %>% 
    select(pft, alldat) %>% 
    unnest()

area_means_wide_all <- area_proc %>% 
    filter(variable == 'mu') %>% 
    select(-variable) %>% 
    mutate(index = factor(index, trait_order)) %>% 
    spread(index, value)

area_means_wide <- select(area_means_wide_all, -pft)

area_glob_mean <- vec2df(colMeans(area_means_wide), pft = 'global', variable = 'mu')
area_glob_cov <- mat2df(cov(area_means_wide), pft = 'global', variable = 'sigma')
area_glob_cor <- mat2df(cor(area_means_wide), diag = FALSE, pft = 'global', variable = 'corr')
area_glob <- bind_rows(area_glob_mean, area_glob_cov, area_glob_cor)
area_all <- full_join(area_glob, area_proc) %>% 
    mutate(pft = factor(pft, levels = c('global', levels(area_proc[['pft']])))) %>% 
    separate(index, c('xvar', 'yvar'), sep = '\\.\\.', remove = FALSE) %>% 
    mutate(xvar = factor(xvar, levels = trait_order),
           yvar = factor(yvar, levels = trait_order))

npft <- n_distinct(area_all$pft)
pfts <- levels(area_all$pft)
pfts <- pfts[pfts %in% area_all$pft]
pft_colors <- c('black', RColorBrewer::brewer.pal(npft - 1, 'Set1'))
names(pft_colors) <- pfts

############################################################
# Bar chart
dir.create('figures')
pdf('figures/area_corrmat.pdf', width = 12, height = 10)
corrmat_plot_dat <- area_all %>% 
    filter(variable == 'corr', !is.na(xvar), !is.na(yvar)) %>% 
    mutate(xvar = factor(xvar, trait_order[-1]),
           yvar = factor(yvar, trait_order[-length(trait_order)]))
ggplot(corrmat_plot_dat) + 
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
area_both <- area_imp2 %>% 
    select(pft, mu_means, Sigma_means) %>% 
    bind_rows(tibble(pft = 'global', 
                     mu_means = list(colMeans(area_means_wide)),
                     Sigma_means = list(cov(area_means_wide))))

ellipsePlot <- function(xvar, yvar, radius = qnorm(0.95)) {
    vrs <- c(xvar, yvar)
    plt <- area_both %>% 
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

trait_combinations <- distinct(area_all %>% filter(variable == 'corr'), xvar, yvar) %>% 
    mutate_all(as.character)

dir.create('figures/pairs_area')
plots <- trait_combinations %$% 
    map2(xvar, yvar, ellipsePlot, radius = qnorm(0.9))
files <- trait_combinations %$%
    file.path('figures/pairs_area', paste(xvar, yvar, 'pdf', sep = '.'))
walk2(files, plots, ggsave, device = pdf, width = 10, height = 8)
