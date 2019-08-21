context("Fit hierarchical model")

test_that("Fit iris data by species", {
  dat_full <- list(
    dat = as.matrix(iris[, -5]),
    groups = factor(iris[["Species"]])
  )

  niter <- 5000
  nchains <- 2
  autofit <- TRUE
  keep_samples <- 2000

  future::plan("multiprocess")
  result <- fit_mvnorm_hier(
    dat_full$dat, dat_full$groups, niter = niter, nchains = nchains,
    autofit = autofit, keep_samples = keep_samples
  )

  summary_proc <- result$summary_table %>%
    dplyr::filter(variable == "Corr") %>%
    dplyr::select(index, group, Mean, `2.5%`, `97.5%`) %>%
    tidyr::separate(index, c("xvar", "yvar"), sep = "\\.\\.") %>%
    tidyr::nest(-xvar, -yvar)

  if (FALSE) {
    bar_plot_matrix(summary_proc, summary_proc, param_order = colnames(iris)[-5], diag_cex = 1)
  }
})
