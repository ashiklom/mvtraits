rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

# Simulate some data
niris <- nrow(iris)
inmat <- cbind(iris[,-5], rnorm(niris), rnorm(niris, 5), rnorm(niris, -3))
colnames(inmat) <- NULL
mu <- colMeans(inmat)
Sig <- cov(inmat)

N <- 5000
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)

# Randomly remove a fraction of the data
dat <- dat_all
nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)
dat[miss] <- NA

dat_filled <- mvtraits:::mvnorm_fill_missing(dat, mu, chol(Sig))
imputed <- dat_filled
imputed[!is.na(dat)] <- NA

print('True mean')
print(mu)
print('Imputed means')
mu_imp <- colMeans(imputed, na.rm = TRUE)
mu_imp
print('Difference')
mu_imp - mu

print('True sigma')
print(Sig)
print('Imputed sigma')
Sig_imp <- cov(imputed, use = 'pairwise.complete.obs')
Sig_imp
print('Difference')
Sig_imp - Sig


if (interactive()) {
    testplot <- function(i, j) {
        plot(dat[,i], dat[,j])
        points(imputed[,i], imputed[,j], col = 'red')
    }
    testplot(2, 4)
}

