rm(list = ls())

library(testthat)
library(mvtraits)

# Fit iris data by species
dat_full <- list(
  dat = as.matrix(iris[, -5]),
  groups = factor(iris[["Species"]])
)

niter <- 5000
nchains <- 2
parallel <- FALSE
autofit <- TRUE
keep_samples <- 2000

message("Running simple hierarchical...")
result <- fit_mvnorm_hier(
  dat_full$dat, dat_full$groups, niter = niter, nchains = nchains,
  parallel = parallel, autofit = autofit, keep_samples = keep_samples
)
message("Done!")

summary_proc <- result$summary_table %>%
  dplyr::filter(variable == "Corr") %>%
  dplyr::select(index, group, Mean, `2.5%`, `97.5%`) %>%
  tidyr::separate(index, c("xvar", "yvar"), sep = "\\.\\.") %>%
  tidyr::nest(-xvar, -yvar)

bar_plot_matrix(summary_proc, summary_proc, param_order = colnames(iris)[-5], diag_cex = 1)
