dat_all <- as.matrix(iris[, -5])
dat <- dat_all
dat[, 4] <- NA
dat[sample(length(dat), 100)] <- NA

dat

N <- 5
seeds <- replicate(N, sample(.Random.seed[-2:0], 7),
                   simplify = FALSE)



future::plan(future.callr::callr())
f <- function(x) {
  s <- sample(1e7, 1)
  mvtraits:::zsetseed(s)
  random_mvnorm(1, c(0, 0), diag(2))
}
furrr::future_map(1:5, f)
## furrr::future_map(seq_len(N), f, .options = furrr::future_options(seed = 43L))
