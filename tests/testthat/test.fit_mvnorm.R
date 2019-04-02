context("Multivariate fit")

dat <- as.matrix(iris[, -5])
niter <- 10000
nchains <- 8
keep_samples <- Inf
max_attempts <- 20
autofit <- TRUE

test_that("Simple multivariate fit, sequentially", {
  result <- fit_mvnorm(
    dat, niter = niter, nchains = nchains,
    autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
  )
})

test_that("Multivariate fit, in parallel", {

  future::plan("multiprocess")
  result <- fit_mvnorm(
    dat, niter = niter, nchains = nchains,
    autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
  )

})

#param <- c("par01", "par02")
#ellipse_dat <- ellipse_axes(
  #mean = result$means$mu[param],
  #cov = result$means$Sigma[param, param]
#)
#draw_ellipse(ellipse_dat)
#draw_majoraxis(ellipse_dat)
