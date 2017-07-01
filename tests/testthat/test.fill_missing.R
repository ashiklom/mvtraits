rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

ss <- function(x, y) {
    sum((y - x)^2)
}

# Simulate some data
rand <- random_data(frac_miss = 0.50)
attach(rand)
nparam <- length(rand$mu)

# Fill in missing data
setup <- setup_missing(dat)
dat_filled <- mvtraits:::c_alt_fill_missing(dat, mu, Sigma, setup)

imputed <- dat_filled
imputed[!is.na(dat)] <- NA

mu_imp <- colMeans(imputed, na.rm = TRUE)
Sigma_imp <- cov(imputed, use = 'pairwise.complete.obs')

test_that('Imputed values are close to true values', {
          expect_less_than(ss(mu_imp, mu), 0.05)
          expect_less_than(ss(mu_imp, mu), 0.01)
})

if (interactive()) {
    testplot <- function(i, j) {
        plot(dat[,i], dat[,j], pch = '.')
        points(imputed[,i], imputed[,j], pch = '.', col = 'red')
        mixtools::ellipse(mu = mu[c(i,j)], sigma = Sigma[c(i,j), c(i,j)], )
    }
    par(mfrow = c(nparam, nparam))
    for (i in seq_len(nparam)) {
        for (j in seq_len(nparam)) {
            if (i == j) {
                plot(0, 0, type = 'n')
            } else {
                testplot(i, j)
            }
        }
    }
}

