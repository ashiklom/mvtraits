#' Fit multivariate normal distribution to data
#'
#' @param dat Input data, as a `data.frame` or matrix
#' @param niter Number of iterations per attempt fitting attempt.
#' @param priors Named list of prior overrides. See [gibbs_default_priors()] for
#'   details.
#' @param inits Named list of initial conditions
#' @param nchains Number of chains (default = 3)
#' @param parallel (Logical) Whether or not to run chains in parallel (default = `TRUE`)
#' @param autofit (Logical) If `TRUE`, keep re-trying until model has converged (default = `FALSE`)
#' @param max_attempts Maximum number of attempts for `autofit` (default = 10)
#' @param keep_samples Maximum number of samples to keep in results (default = Inf)
#' @param threshold Gelman PSRF threshold for convergence (default = 1.15)
#' @param save_progress If not `NULL`, path to file in which to save current
#' @param progress (Logical) If `TRUE`, display a progress bar while sampling
#'   chain status between fitting attempts. If `NULL` (default), default to
#'   `TRUE` if using a sequential `future::plan` and `FALSE`
#' @return List containing a tidy summary table of the fitted parameters
#'   (`summary_table`), a list of the summary statistics (`stats`), and the raw
#'   MCMC samples as an `mcmc.list` (`samples`).
#' @useDynLib mvtraits
#' @export
fit_mvnorm <- function(dat, niter = 5000, priors = list(), inits = list(), nchains = 3,
                       autofit = FALSE, max_attempts = 10, keep_samples = Inf,
                       threshold = 1.15, save_progress = NULL,
                       progress = NULL) {

  # Automatically determine whether or not to print a progress bar
  if (is.null(progress)) {
    progress <- inherits(future::plan(), "sequential")
  }

  chainseq <- seq_len(nchains)

  nparam <- ncol(dat)
  param_names <- colnames(dat)
  if (is.null(param_names)) {
    param_names <- sprintf("par%02d", seq_len(nparam))
  }

  ndat <- nrow(dat)
  setup <- setup_missing(dat)

  # Where missing, use default priors
  default_priors <- gibbs_default_priors(nparam)
  if (!is.null(priors)) {
    priors <- modifyList(default_priors, priors)
  } else {
    priors <- default_priors
  }

  # Set priors in environment
  mu0 <- priors[["mu_global"]]
  Sigma0 <- priors[["Sigma_global"]]
  v0 <- priors[["v_global"]]
  S0 <- priors[["S_global"]]

  Sigma0_inv <- solve(Sigma0)

  # Draw default initial conditions from priors
  mu <- list()
  Sigma <- list()
  for (n in chainseq) {
    mu[[n]] <- random_mvnorm(1, mu0, Sigma0)[1, ]
    names(mu[[n]]) <- param_names
    Sigma[[n]] <- solve(rWishart(1, v0 + nparam + 1, S0)[, , 1])
    dimnames(Sigma[[n]]) <- list(param_names, param_names)
  }
  default_inits <- list(mu = mu, Sigma = Sigma)
  if (!is.null(inits)) {
    inits <- modifyList(default_inits, inits)
  } else {
    inits <- default_inits
  }

  sampler <- list(
    fun = sample_mvnorm,
    init_fun = function(n, inits) {
      list(
        mu = inits[["mu"]][[n]],
        Sigma = inits[["Sigma"]][[n]]
      )
    },
    args = list(
      niter = niter,
      dat = dat,
      mu0 = mu0,
      Sigma0_inv = Sigma0_inv,
      v0 = v0,
      S0 = S0,
      setup = setup,
      progress = progress
    )
  )

    message("Running sampler...")
    raw_samples <- run_until_converged(
      sampler = sampler,
      model_type = 'multi',
      inits = inits,
      nchains = nchains,
      max_attempts = max_attempts,
      save_progress = save_progress,
      threshold = threshold,
      keep_samples = keep_samples,
      autofit = autofit
    )

  message("Calculating correlation matrices...")
  raw_samples_corr <- add_correlations(raw_samples)
  rm(raw_samples)

  message("Converting samples to coda mcmc.list object...")
  samples_mcmc <- results2mcmclist(raw_samples_corr, type = "multi")
  rm(raw_samples_corr)

  niter <- coda::niter(x = samples_mcmc)

  message("Preparing summary table...")
  summary_table <- summary_df(
    window(samples_mcmc, start = floor(niter / 2)),
    group = NULL
  )

  stats <- c("Mean", "2.5%", "97.5%")
  mu_stats <- sapply(
    stats,
    function(x) summary2vec(summary_table, x, variable == "mu"),
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  Sigma_stats <- sapply(
    stats,
    function(x) summary2mat(summary_table, x, variable == "Sigma"),
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  Corr_stats <- sapply(
    stats,
    function(x) summary2mat(summary_table, x, variable == "Corr"),
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  list(
    summary_table = summary_table,
    stats = list(
      mu = mu_stats,
      Sigma = Sigma_stats,
      Corr = Corr_stats
    ),
    samples = samples_mcmc
  )
}
