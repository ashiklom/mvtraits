context("Bootstrap missing data")

dat <- as.matrix(iris[, -5])
groups <- iris[, 5]
set.seed(1234123)

# Make 25% of the data go missing
ndata <- prod(dim(dat))
nmissing <- floor(ndata * 0.25)
imissing <- sample(ndata, nmissing)
y <- dat
ytrue <- y[imissing]
y[imissing] <- NA

cols <- c("Sepal.Length", "Sepal.Width",
          "Petal.Length", "Petal.Width")

test_that("Impute with multivariate normal fit", {
  result <- fit_mvnorm(y)
  out <- bootstrap_missing(result, y)
  expect_equal(colnames(out), cols)
  means <- apply(out, c(1, 2), mean)[imissing]
  lo <- apply(out, c(1, 2), quantile, probs = 0.025)
  hi <- apply(out, c(1, 2), quantile, probs = 0.975)
  rmse <- sqrt(mean((ytrue - means)^2))
  # RMSE should pretty consistently be below 0.6
  expect_lt(rmse, 0.6)
  # 1:1 regression should be pretty close to 1
  fit <- lm(ytrue ~ means)
  fitcoef <- coefficients(fit)
  expect_lt(abs(fitcoef[1]), 0.1)
  expect_lt(abs(fitcoef[2] - 1), 0.05)
})

test_that("Impute with hierarchical fit", {
  result <- fit_mvnorm_hier(y, groups)
  out <- bootstrap_missing_hier(result, y, groups = groups)
  expect_equal(colnames(out), cols)
  means <- apply(out, c(1, 2), mean)[imissing]
  lo <- apply(out, c(1, 2), quantile, probs = 0.025)
  hi <- apply(out, c(1, 2), quantile, probs = 0.975)
  rmse <- sqrt(mean((ytrue - means) ^ 2))
  expect_lt(rmse, 0.6)
  fit <- lm(ytrue ~ means)
  fitcoef <- coefficients(fit)
  expect_lt(abs(fitcoef[1]), 0.1)
  expect_lt(abs(fitcoef[2] - 1), 0.05)
})
