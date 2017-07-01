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
Sigma <- cov(inmat)

N <- 5000
dat_all <- mvtraits:::c_random_mvnorm(N, t(mu), Sigma)

# Randomly remove a fraction of the data
dat <- dat_all
nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)
dat[miss] <- NA
#scramble <- 4:1
#scramble <- c(4,1,2,3)
#dat <- dat[,scramble]
#mu_scr <- mu[scramble]
#Sigma_scr <- Sigma[scramble, scramble]

setup <- setup_missing(dat)
dat_filled <- mvtraits:::c_alt_fill_missing(dat, mu, Sigma, setup)

imputed <- dat_filled
imputed[!is.na(dat)] <- NA

print('True mean')
print(mu)
print('Imputed means')
mu_imp <- colMeans(imputed, na.rm = TRUE)
print(mu_imp)
print('Difference')
print(mu_imp - mu)

print('True sigma')
print(Sigma)
print('Imputed sigma')
Sigma_imp <- cov(imputed, use = 'pairwise.complete.obs')
print(Sigma_imp)
print('Difference')
print(Sigma_imp - Sigma)

if (interactive()) {
    testplot <- function(i, j) {
        plot(dat[,i], dat[,j])
        points(imputed[,i], imputed[,j], col = 'red')
        #inds <- scramble[c(i, j)]
        inds <- c(i,j)
        mixtools::ellipse(mu = mu[inds], sigma = Sigma[inds, inds], )
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

