context("Bootstrap missing data")

dat <- as.matrix(iris[, -5])
groups <- iris[, 5]
set.seed(1234123)

# Make 25% of the data go missing
ndata <- prod(dim(dat))
nmissing <- floor(ndata * 0.25)
imissing <- sort(sample(ndata, nmissing))
y <- dat
ytrue <- y[imissing]
y[imissing] <- NA

cols <- c("Sepal.Length", "Sepal.Width",
          "Petal.Length", "Petal.Width")

test_that("Impute with multivariate normal fit", {
  result <- fit_mvnorm(y)
  out <- bootstrap_missing(result, y)
  expect_equal(colnames(out), cols)
  expect_equal(attr(out, "imissing"), which(is.na(y), arr.ind = TRUE))
  out_smry <- tidy_bootstrap(out)
  means <- out_smry$Mean
  rmse <- sqrt(mean((ytrue - means)^2))
  expect_lt(rmse, 0.6)
  fit <- lm(ytrue ~ means)
  fitcoef <- coefficients(fit)
  expect_lt(abs(fitcoef[1]), 0.1)
  expect_lt(abs(fitcoef[2] - 1), 0.05)
})

test_that("Impute with hierarchical fit", {
  result <- fit_mvnorm_hier(y, groups)
  out <- bootstrap_missing_hier(result, y, groups = groups)
  expect_equal(colnames(out), cols)
  expect_equal(attr(out, "imissing"), which(is.na(y), arr.ind = TRUE))
  out_smry <- tidy_bootstrap(out)
  rmse <- sqrt(mean((ytrue - out_smry$Mean) ^ 2))
  expect_lt(rmse, 0.6)
  fit <- lm(ytrue ~ out_smry$Mean)
  fitcoef <- coefficients(fit)
  expect_lt(abs(fitcoef[1]), 0.1)
  expect_lt(abs(fitcoef[2] - 1), 0.05)
})
