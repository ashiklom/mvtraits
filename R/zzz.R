.onLoad <- function(libname, pkgname) {
  z <- runif(1)   # Necessary to initialize .Random.seed variable
  s <- .Random.seed[-2:0]
  seed <- sample(s, 1)
  RcppZiggurat::zsetseed(seed)
}
