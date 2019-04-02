context("Matrix helpers")

test_that("Flatten-unflatten matrix is reversible", {
  cv <- cov(iris[, -5])
  flat <- flatten_matrix(cv)
  unflat <- lowerdiag2mat(flat)

  expect_equal(cv, unflat)
})
