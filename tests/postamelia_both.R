library(Amelia)
library(tidyverse)
library(magrittr)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

area_imp <- readRDS('data/amelia_pft.rds')

trait_order_area <- c('leaf_lifespan', 'LMA', 'Narea', 'Parea', 'Rdarea', 'Vcmax_area', 'Jmax_area')
area_imp2 <- imp2(area_imp, trait_order_area)
area_proc <- proc(area_imp2)
area_means_wide <- means_wide(area_proc, trait_order_area)
area_all <- bind_all(area_means_wide, area_proc, trait_order_area)

mass_imp <- readRDS('data/amelia_pft_mass.rds')

trait_order_mass <- c('leaf_lifespan', 'SLA', 'Nmass', 'Pmass', 'Rdmass', 'Vcmax_mass', 'Jmax_mass')
mass_imp2 <- imp2(mass_imp, trait_order_mass)
mass_proc <- proc(mass_imp2)
mass_means_wide <- means_wide(mass_proc, trait_order_mass)
mass_all <- bind_all(mass_means_wide, mass_proc, trait_order_mass)

npft <- n_distinct(area_all$pft)
pfts <- levels(area_all$pft)
pfts <- pfts[pfts %in% area_all$pft]
pft_colors <- c('black', RColorBrewer::brewer.pal(npft - 1, 'Set1'))
names(pft_colors) <- pfts

both_all <- bind_rows(area_all, mass_all)

############################################################
# Bar chart
dir.create('figures')
pdf('figures/both_corrmat.pdf', width = 12, height = 10)
corrmat_plot_dat <- area_all %>% 
    filter(variable == 'corr', !is.na(xvar), !is.na(yvar)) %>% 
    mutate(xvar = factor(xvar, trait_order_area[-1]),
           yvar = factor(yvar, trait_order_area[-length(trait_order_area)]))
ggplot(corrmat_plot_dat) + 
    aes(x = pft, y = value, fill = pft) + 
    geom_col() + 
    facet_wrap(~xvar + yvar, drop = FALSE) + 
    scale_fill_manual(values = pft_colors) + 
    theme(axis.text.x = element_blank()) + 
    xlab('Plant functional type') + 
    ylab('Mean correlation coefficient')

corrmat_plot_dat <- mass_all %>% 
    filter(variable == 'corr', !is.na(xvar), !is.na(yvar)) %>% 
    mutate(xvar = factor(xvar, trait_order_mass[-1]),
           yvar = factor(yvar, trait_order_mass[-length(trait_order_mass)]))
ggplot(corrmat_plot_dat) + 
    aes(x = pft, y = value, fill = pft) + 
    geom_col() + 
    facet_wrap(~yvar + xvar, drop = FALSE) + 
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
