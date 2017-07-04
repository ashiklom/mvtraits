library(tidyverse)

raw_mass <- readRDS('summaries/mass_hier.rds')
raw_area <- readRDS('summaries/area_hier.rds')

try_data <- readRDS('~/Projects/new-phytologist-traits/preprocess-try/traits_analysis.rds')
pft_levels <- try_data %>% pull(pft) %>% levels() %>% c('global', .)
group_names <- raw_mass %>% distinct(group) %>% pull(group) %>% sort()

proc <- . %>%
    mutate(pft = factor(group, levels = group_names) %>% 'levels<-'(pft_levels)) %>% 
    separate(index, c('xvar', 'yvar'), sep = '\\.\\.', extra = 'merge')

mass_proc <- raw_mass %>% proc
area_proc <- raw_area %>% proc %>% 
    mutate(yvar = gsub('group..\\.\\.', '', yvar))

# Assign PFT colors
pft_colors <- c('black', RColorBrewer::brewer.pal(length(pft_levels) - 1, 'Paired'))
names(pft_colors) <- pft_levels

grid_corr_plot <- function(dat) {
    ggplot(dat %>% filter(variable == 'Corr')) +
        aes(x = pft, y = Mean, ymin = `2.5%`, ymax = `97.5%`, fill = pft) + 
        geom_col() +
        geom_errorbar(width = 0.2) +
        facet_wrap(~xvar + yvar, drop = FALSE) + 
        scale_fill_manual(values = pft_colors) + 
        theme(axis.text.x = element_blank())
}

fig_dir <- file.path('figures', 'gibbs')
dir.create(fig_dir, showWarnings = FALSE)

mass_plt <- grid_corr_plot(mass_proc)
ggsave(filename = file.path(fig_dir, 'mass_corrmat.pdf'), plot = mass_plt, width = 10, height = 8)
area_plt <- grid_corr_plot(area_proc)
ggsave(filename = file.path(fig_dir, 'area_corrmat.pdf'), plot = area_plt, width = 10, height = 8)

mass_proc %>% 
    filter(pft == 'global', variable == 'Sigma') %>% 
    ggplot() + 
        aes(x = 1, y = Mean, ymin = `2.5%`, ymax = `97.5%`) + 
        geom_col() + 
        geom_errorbar(width = 0.2) +
        facet_wrap(~xvar + yvar, drop = FALSE)
