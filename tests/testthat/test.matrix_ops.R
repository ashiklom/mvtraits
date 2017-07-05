library(testthat)
library(mvtraits)

cv <- cov(iris[,-5])
flat <- flatten_matrix(cv)
unflat <- lowerdiag2mat(flat)

test_that('Flatten-unflatten matrix is reversible', {
          expect_equal(cv, unflat)
})
