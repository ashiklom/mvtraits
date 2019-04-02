context("Fill missing data")

ss <- function(x, y) {
  sum((y - x) ^ 2)
}

# Simulate some data
rand <- random_data_multi(frac_miss = 0.50)
nparam <- length(rand$mu)

# Fill in missing data
dat_filled <- mvnorm_fill_missing(rand$dat, rand$mu, rand$sigma)

imputed <- dat_filled
imputed[!is.na(rand$dat)] <- NA

mu_imp <- colMeans(imputed, na.rm = TRUE)
Sigma_imp <- cov(imputed, use = 'pairwise.complete.obs')

test_that("Imputed values are close to true values", {
  expect_lt(ss(mu_imp, rand[["mu"]]), 0.05)
  expect_lt(ss(Sigma_imp, rand[["sigma"]]), 0.1)
})

if (exists('doplot')) {
    testplot <- function(i, j) {
        plot(dat[,i], dat[,j], pch = '.')
        points(imputed[,i], imputed[,j], pch = '.', col = 'red')
        mixtools::ellipse(mu = mu[c(i,j)], sigma = sigma[c(i,j), c(i,j)], )
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

