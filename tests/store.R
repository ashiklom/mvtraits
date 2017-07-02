#library(mvtraits)
devtools::load_all('.')
library(testthat)

m <- cov(iris[,-5])

mb <- microbenchmark::microbenchmark(store_covmat(m), store_covmat2(m), times = 10000)
print(mb)
