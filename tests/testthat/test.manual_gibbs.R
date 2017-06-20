# Data
devtools::load_all('.')

set.seed(666)
nx <- 10000
nmiss <- 2000
irismat <- as.matrix(iris[,-5])

true_mu <- colMeans(irismat)
true_Sigma <- cov(irismat)
nparam <- length(true_mu)
x <- mvtnorm::rmvnorm(nx, true_mu, true_Sigma)
xobs <- x
xobs[sample.int(length(xobs), nmiss)] <- NA

samples <- fit_mvnorm(dat = xobs)

print('Mu:')
mu_mle <- colMeans(samples$mu)
mu_mle
apply(samples$mu, 2, sd)

print("True mu:")
true_mu

print("Mu diff:")
mu_mle - true_mu

print('Sigma:')
sigma_mle <- apply(samples$Sigma, 2:3, mean)
sigma_mle
apply(samples$Sigma, 2:3, sd)

omega_samples <- array(NA_real_, dim(samples$Sigma))
for (i in seq_len(dim(samples$Sigma)[1])){
    omega_samples[i,,] <- cov2cor(samples$Sigma[i,,])
}
omega_mle <- apply(omega_samples, 2:3, mean)

print('True sigma')
true_Sigma

print("Sigma diff:")
sigma_mle - true_Sigma
