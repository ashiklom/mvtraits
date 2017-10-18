library(testthat)
library(mvtraits)

rand <- random_data_multi()

niter <- 10000
nchains <- 2
parallel <- FALSE
keep_samples <- Inf
max_attempts <- 20
autofit <- TRUE

message("Running simple multivariate...")
result <- fit_mvnorm(
  rand$dat, niter = niter, nchains = nchains, parallel = parallel,
  autofit = autofit, keep_samples = keep_samples, max_attempts = max_attempts
)
message("Done!")
# R version: 25 seconds
# R main, C functions: 11 seconds
# All C:

param <- c("par01", "par02")
ellipse_dat <- ellipse_axes(
  mean = result$means$mu[param],
  cov = result$means$Sigma[param, param]
)
draw_ellipse(ellipse_dat)
draw_majoraxis(ellipse_dat)
