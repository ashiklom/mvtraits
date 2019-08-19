context("Fit multivariate model")

dat <- as.matrix(iris[,-5])

niter <- 10000
nchains <- 2
keep_samples <- Inf
max_attempts <- 20
autofit <- TRUE

test_that("Fit multivariate normal sequentially", {
  result <- fit_mvnorm(
    dat, niter = niter, nchains = nchains, parallel = FALSE,
    autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
  )
})

test_that("Fit multivariate normal in parallel", {
  result <- fit_mvnorm(
    dat, niter = niter, nchains = nchains, parallel = TRUE,
    autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
  )
})
