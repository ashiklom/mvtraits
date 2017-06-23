rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

# Simulate some data
niris <- nrow(iris)
#inmat <- cbind(iris[,-5], rnorm(niris), rnorm(niris, 5), rnorm(niris, -3))
#colnames(inmat) <- NULL
inmat <- as.matrix(iris[,-5])
nparam <- ncol(inmat)
mu <- colMeans(inmat)
Sig <- cov(inmat)

N <- 5000
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)

# Randomly remove a fraction of the data
dat <- dat_all
nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)
dat[miss] <- NA
scramble <- 4:1
#scramble <- c(4,1,2,3)
dat <- dat[,scramble]
mu_scr <- mu[scramble]
Sig_scr <- Sig[scramble, scramble]

dat_filled <- mvtraits:::mvnorm_fill_missing(dat, mu_scr, chol(Sig_scr))
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
        inds <- scramble[c(i, j)]
        mixtools::ellipse(mu = mu[inds], sigma = Sig[inds, inds], )
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

