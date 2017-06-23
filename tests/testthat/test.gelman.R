ignore <- function() {
mu <- colMeans(iris[,-5])
Sig <- cov(iris[,-5])
Cor <- cov2cor(Sig)

N <- 500
ngroup <- 7
dat_all <- mvtnorm::rmvnorm(N, mu, Sig)
groups <- sample.int(ngroup, N, replace = TRUE)

# Randomly remove a fraction of the data
dat <- dat_all

#nmiss <- round(length(dat) * 0.5)
#miss <- sample.int(length(dat), size = nmiss)

#dat[miss] <- NA

# Constants
nparam <- ncol(dat)
n <- nrow(dat)

# Set priors
mu0 <- rep(0, nparam)
k0 <- 1
v0 <- 0
S0 <- diag(1, nparam)
Sinv0 <- solve(S0)

# Constant parameters
k_n <- k0 + n
v_n <- v0 + n

# Update parameters
xbar <- colMeans(dat)
mu_n <- (n / k_n) * xbar + (k0 / k_n) * mu0
Q0 <- cov(dat) * n
xbm0 <- (xbar - mu0) %*% t(xbar - mu0)
Sinv_n <- Sinv0 + Q0 + (k0 * n) / k_n * xbm0
S_n <- solve(Sinv_n)

# Draw values
Sigma_inv <- rWishart(1, v_n + nparam + 1, S_n)[,,1]
Sigma <- solve(Sigma_inv)
mu <- mvtnorm::rmvnorm(1, mu_n, Sigma / k_n)
}
