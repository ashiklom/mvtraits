library(testthat)
library(mvtraits)
library(mvtnorm)

# Simulate some data
mu <- c(10, 5, -5)
Sig <- matrix(c(1, 0.7, 0.8,
                0.7, 2, 0.2,
                0.8, 0.2, 1.5), 
              3, 3)
N <- 100
dat <- rmvnorm(N, mu, Sig)

# Randomly remove some data
nmiss <- 20
miss <- sample.int(length(dat), size = nmiss)
dat[miss] <- NA
dat[,3] <- NA
dat <- cbind(dat, "pft" = sample.int(5, N, replace = TRUE))
dat <- as.data.table(dat)
dir.create("output", showWarnings = FALSE)

fit_uni <- runModel("uni", dat, NA)
saveRDS(fit_uni, "output/uni.rds")
fit_uni_1 <- runModel("uni", dat, 1)
saveRDS(fit_uni_1, "output/uni.1.rds")
fit_multi <- runModel("multi", dat, NA)
saveRDS(fit_multi, "output/multi.rds")
fit_multi_1 <- runModel("multi", dat, 1)
saveRDS(fit_multi_1, "output/multi.1.rds")
fit_hier <- runModel("hier", dat, NA)
saveRDS(fit_hier, "output/hier.rds")
