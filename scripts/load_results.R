arg <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arg) == 1)

fname <- arg[1]

library(mvtraits)
library(tidyverse)

result_raw <- readRDS(fname)

col_names <- colnames(result_raw[[1]][['mu_group']])
grps_split <- strsplit(col_names, '\\.\\.')
grps <- sapply(grps_split, '[', 1)
ugroups <- unique(grps)
ngroups <- length(ugroups)
niter <- nrow(result_raw[[1]][['mu_global']])

burn <- function(chain, start = floor(niter / 2)) {
    lapply(chain, function(x) x[-start:0,])
}

result_burned <- lapply(result_raw, burn)
result_corr <- add_correlations(result_burned, hier = TRUE, ngroups = ngroups)
result_mcmc <- results2mcmclist(result_corr, 'hier')
result_summary <- summary_df(result_mcmc, group = TRUE)

summary_dir <- 'summaries'
dir.create(summary_dir, showWarnings = FALSE)
save_fname <- file.path(summary_dir, basename(fname))
saveRDS(result_summary, save_fname)
