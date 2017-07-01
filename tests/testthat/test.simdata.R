rm(list = ls())

library(testthat)
if (interactive()) {
    devtools::load_all('.')
} else {
    library(mvtraits)
}

# Simulate some data
rand <- random_data()
attach(rand)

custom_inputs <- list()

#fit_uni <- runModel('uni', dat[groups == 1,], iter = 100, max.attempts = 1)
#fit_multi <- runModel('multi', dat[groups == 1,], iter = 100, max.attempts = 1)
#fit_hier <- runModel('hier', dat, groups = groups, iter = 100, max.attempts = 1)

#file.remove(list.files('testmodel_*.rds'))
