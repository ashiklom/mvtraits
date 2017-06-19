# Data
source('R/manual_gibbs.R')

set.seed(666)
nparam <- 4
nx <- 10000
nmiss <- 30000
x <- cbind(rnorm(nx, -4), rnorm(nx, -2), rnorm(nx, 2), rnorm(nx, 4))
xobs <- x
xobs[sample.int(length(xobs), nmiss)] <- NA

# Prior parameters
mu0 <- rep(0, nparam)
Sigma_0 <- diag(10, nparam)
Sigma_0_inv <- solve(Sigma_0)
v0 <- nparam
S0 <- diag(nparam)

# Initial conditions
mu <- rep(0, nparam)
Sigma <- diag(nparam)

ngibbs <- 1000
mu_samp <- matrix(NA_real_, ngibbs, nparam)
Sigma_samp <- array(NA_real_, c(ngibbs, nparam, nparam))

# Begin sampler
Rprof()
pb <- txtProgressBar(1, ngibbs, style = 3)
for (i in seq_len(ngibbs)) {
    setTxtProgressBar(pb, i)
    Sigma_inv <- solve(Sigma)
    Sigma_chol <- chol(Sigma)
    y <- mvnorm_fill_missing(xobs, mu, Sigma_chol)
    ybar <- colMeans(y)
    mu <- draw_mu(ybar, nx, Sigma_inv, mu0, Sigma_0_inv)
    Sigma <- draw_Sigma(x, mu, v0, S0)
    # Store outputs
    mu_samp[i,] <- mu
    Sigma_samp[i,,] <- Sigma
}
close(pb)
Rprof(NULL)

print('Mu:')
colMeans(mu_samp)
apply(mu_samp, 2, sd)

print('Sigma:')
apply(Sigma_samp, 2:3, mean)
apply(Sigma_samp, 2:3, sd)
