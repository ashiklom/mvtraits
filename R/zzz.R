.onLoad <- function(libname, pkgname) {
  z <- runif(1)   # Necessary to initialize .Random.seed variable
  zsetseed(sample(1e7, 1))
}
