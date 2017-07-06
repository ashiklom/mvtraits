library(tidyverse)
library(mvtraits)
if (interactive()) devtools::install()

raw_mass <- readRDS('summaries/mass_hier.rds')
raw_area <- readRDS('summaries/area_hier.rds')

try_data <- readRDS('~/Projects/new-phytologist-traits/preprocess-try/traits_analysis.rds')
pft_levels <- try_data %>% pull(pft) %>% levels() %>% c('global', .)
group_names <- raw_mass %>% distinct(group) %>% pull(group) %>% sort()

proc <- function(dat, params) {
    dat %>% 
        select(-rowname) %>% 
        mutate(index = gsub('group..\\.\\.', '', index)) %>% 
        rearrange_df(params) %>% 
        mutate(pft = factor(group, group_names) %>% 'levels<-'(pft_levels))
}

mass_params <- c('leaf_lifespan', 'SLA', 'Nmass', 'Pmass', 'Rdmass', 'Vcmax_mass', 'Jmax_mass')
area_params <- gsub('mass$', 'area', mass_params)

proc_mass <- proc(raw_mass, mass_params)
proc_area <- proc(raw_area, area_params)

# Assign PFT colors
pft_colors <- c('black', RColorBrewer::brewer.pal(length(pft_levels) - 1, 'Paired'))
names(pft_colors) <- pft_levels

grid_corr_plot <- function(dat, ...) {
    summary_barplot(dat, 'Corr', ...) + scale_fill_manual(values = pft_colors)
}

fig_dir <- file.path('figures', 'gibbs')
dir.create(fig_dir, showWarnings = FALSE)

mass_plt <- grid_corr_plot(proc_mass)
ggsave(filename = file.path(fig_dir, 'mass_corrmat.pdf'), plot = mass_plt, width = 10, height = 8)

area_plt <- grid_corr_plot(proc_area, facet_list = list(facets = ~yparam + xparam))
ggsave(filename = file.path(fig_dir, 'area_corrmat.pdf'), plot = area_plt, width = 10, height = 8)
