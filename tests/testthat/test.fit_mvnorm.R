library(testthat)
library(mvtraits)

dat <- as.matrix(iris[,-5])

niter <- 10000
nchains <- 2
parallel <- FALSE
keep_samples <- Inf
max_attempts <- 20
autofit <- TRUE

message("Running simple multivariate...")
result <- fit_mvnorm(
  dat, niter = niter, nchains = nchains, parallel = parallel,
  autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
)
message("Done!")

message("Running simple multivariate in parallel...")
parallel <- TRUE
result <- fit_mvnorm(
  dat, niter = niter, nchains = nchains, parallel = parallel,
  autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
)
message("Done!")
