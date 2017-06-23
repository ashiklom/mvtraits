nsamp <- 5000
nparam <- 5 

target_mean <- 1
target_sd <- 1

v0 <- 20
df_n <- v0 + nparam + 1

diag_s <- target_mean / df_n
S0 <- diag(diag_s, nparam)
samp <- rWishart(nsamp, df_n, S0)

x <- Wishart_prior_param(1, 0.01, 5)
samp <- rWishart(nsamp, x$v0 + nparam + 1, x$S0)

Mean <- apply(samp, 1:2, mean)
Mean
Var <- apply(samp, 1:2, sd)
Var

invsamp <- samp
for (i in seq_len(nsamp)) {
    invsamp[,,i] <- solve(samp[,,i])
}
invMean <- apply(invsamp, 1:2, mean)
invSD <- apply(invsamp, 1:2, sd)
