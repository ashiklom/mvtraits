rm(list = ls())

library(testthat)
library(mvtraits)
#devtools::load_all('.')
library(mvtnorm)
library(clusterGeneration)

# Simulate some data
mu <- c(10, 5, 0, -5, 10)
Sig <- genPositiveDefMat(length(mu))$Sigma

N <- 1000
ngroup <- 7
dat_all <- rmvnorm(N, mu, Sig)
groups <- sample.int(ngroup, N, replace = TRUE)

# Randomly remove some data
#dat <- dat_all[groups == 1,]
dat <- dat_all

nmiss <- round(length(dat) * 0.5)
miss <- sample.int(length(dat), size = nmiss)

dat[miss] <- NA

custom_inputs <- list()

#dat[dat[,"pft"] == 1, "pft"] <- 2
#dat[1, "pft"] <- 1
#dat <- as.data.table(dat)
dir.create("output", showWarnings = FALSE)

fit_uni <- runModel('uni', dat[groups == 1,])
fit_multi <- runModel('multi', dat[groups == 1,])
fit_hier <- runModel('hier', dat, groups = groups, iter = 100)
#fit_uni <- runModel("uni", dat, NA)
#saveRDS(fit_uni, "output/uni.rds")
#fit_uni_1 <- runModel("uni", dat, 1)
#saveRDS(fit_uni_1, "output/uni.1.rds")
#fit_multi <- runModel("multi", dat, NA)
#saveRDS(fit_multi, "output/multi.rds")
#for (i in seq_len(npft)) {
    #message("Running pft ", i)
    #fit_multi_pft <- runModel("multi", dat, i)
#}
#saveRDS(fit_multi_1, "output/multi.1.rds")
#fit_hier <- runModel("hier", dat, NA)
#saveRDS(fit_hier, "output/hier.rds")
