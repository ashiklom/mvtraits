library(testthat)
library(mvtraits)
library(mvtnorm)
library(clusterGeneration)

# Simulate some data
mu <- c(10, 5, 0, -5, 10)
Sig <- genPositiveDefMat(length(mu))$Sigma

N <- 1000
npft <- 7
dat <- rmvnorm(N, mu, Sig)

# Randomly remove some data
nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)
dat[miss] <- NA
dat <- cbind(dat, "pft" = sample.int(npft, N, replace = TRUE))
dat[dat[,"pft"] == 1, 3] <- NA
#dat[dat[,"pft"] == 1, "pft"] <- 2
#dat[1, "pft"] <- 1
#dat <- as.data.table(dat)
dir.create("output", showWarnings = FALSE)

#fit_uni <- runModel("uni", dat, NA)
#saveRDS(fit_uni, "output/uni.rds")
#fit_uni_1 <- runModel("uni", dat, 1)
#saveRDS(fit_uni_1, "output/uni.1.rds")
fit_multi <- runModel("multi", dat, NA)
#saveRDS(fit_multi, "output/multi.rds")
#for (i in seq_len(npft)) {
    #message("Running pft ", i)
    #fit_multi_pft <- runModel("multi", dat, i)
#}
#saveRDS(fit_multi_1, "output/multi.1.rds")
fit_hier <- runModel("hier", dat, NA)
#saveRDS(fit_hier, "output/hier.rds")
