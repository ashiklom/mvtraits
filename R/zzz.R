.onLoad <- function(libname, pkgname) {
    z <- runif(1)   # Necessary to initialize .Random.seed variable
    set_R_seed()
}
