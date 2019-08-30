#' Fit a hierarchical multivariate model to data
#'
#' @param groups A character, integer, or factor of group labels, with length
#'   `nrow(dat)`.
#' @inheritParams fit_mvnorm
#' @inherit fit_mvnorm return
#' @export
fit_mvnorm_hier <- function(dat,
                            groups,
                            niter = 5000,
                            priors = list(),
                            inits = list(),
                            nchains = 3,
                            autofit = FALSE,
                            max_attempts = 10,
                            keep_samples = Inf,
                            threshold = 1.15,
                            save_progress = NULL,
                            progress = NULL) {

  if (is.null(progress)) {
    progress <- inherits(future::plan(), "sequential")
  }

  stopifnot(is.matrix(dat), length(groups) == nrow(dat))

  chainseq <- seq_len(nchains)

  nparam <- ncol(dat)
  param_names <- colnames(dat)
  if (is.null(param_names)) {
    param_names <- sprintf("par%02d", seq_len(nparam))
  }

  ngroup <- length(unique(groups))
  if (is.character(groups)) {
    groups <- factor(groups)
  }
  if (is.factor(groups)) {
    group_names <- levels(groups)
  } else {
    group_names <- sprintf("group%02d", seq_len(ngroup))
  }
  igroups <- as.integer(groups)
  ugroups <- sort(unique(igroups))

  setup_bygroup <- lapply(
    ugroups,
    function(x) setup_missing(dat[igroups == x, ])
  )

  # Where missing, use default priors
  default_priors <- gibbs_default_priors(nparam, ngroup)
  if (!is.null(priors)) {
    priors <- modifyList(default_priors, priors)
    if (length(priors) != length(default_priors)) {
      stop(
        "Length of priors (", length(priors), ") ",
        "does not equal length of default priors (",
        length(default_priors), "). ",
        "There is likely a typo in your `prior` name.\n",
        "names(priors): ", paste(names(priors), collapse = ", "),
        "\n",
        "names(default_priors): ", paste(names(default_priors), collapse = ", ")
      )
    }
  } else {
    priors <- default_priors
  }

  # Set priors in environment
  mu0_global <- priors[["mu_global"]]
  Sigma0_global <- priors[["Sigma_global"]]
  v0_global <- priors[["v_global"]]
  S0_global <- priors[["S_global"]]

  mu0_group <- priors[["mu_group"]]
  Sigma0_group <- priors[["Sigma_group"]]
  v0_group <- priors[["v_group"]]
  S0_group <- priors[["S_group"]]

  # Precalculate certain quantities
  Sigma0_global_inv <- solve(Sigma0_global)
  Sigma0_group_inv <- vapply(
    seq_len(ngroup),
    function(x) solve(Sigma0_group[x,,]),
    Sigma0_global_inv
  )
  Sigma0_group_inv <- aperm(Sigma0_group_inv, c(3, 1, 2))

  # Draw initial conditions from priors
  mu_global <- list()
  Sigma_global <- list()
  mu_group <- list()
  Sigma_group <- list()
  for (n in chainseq) {
    mu_global[[n]] <- random_mvnorm(1, mu0_global, Sigma0_global)[1, ]
    names(mu_global[[n]]) <- param_names
    Sigma_global[[n]] <- solve(rWishart(1, v0_global + nparam + 1,
                                        S0_global)[,,1])
    dimnames(Sigma_global[[n]]) <- list(param_names, param_names)

    mu_group[[n]] <- matrix(NA_real_, nrow = ngroup, ncol = nparam)
    dimnames(mu_group[[n]]) <- list(group_names, param_names)
    Sigma_group[[n]] <- array(NA_real_, c(ngroup, nparam, nparam))
    dimnames(Sigma_group[[n]]) <- list(group_names, param_names, param_names)
    for (i in seq_len(ngroup)) {
      mu_group[[n]][i, ] <- random_mvnorm(1, mu0_group[i,], Sigma0_group[i,,])
      #Sigma_group[[n]][i,,] <- solve(rWishart(1, v0_group[i] + nparam + 1, S0_group[i,,])[,,1])
      Sigma_group[[n]][i, , ] <- diag(1, nparam)
    }
  }
  default_inits <- list(mu_global = mu_global,
                        Sigma_global = Sigma_global,
                        mu_group = mu_group,
                        Sigma_group = Sigma_group)
  if (!is.null(inits)) {
    inits <- modifyList(default_inits, inits)
  } else {
    inits <- default_inits
  }

  sampler <- list(
    fun = sample_mvnorm_hier,
    init_fun = function(n, inits) {
      list(
        mu_global = inits[["mu_global"]][[n]],
        Sigma_global = inits[["Sigma_global"]][[n]],
        mu_group = inits[["mu_group"]][[n]],
        Sigma_group = inits[["Sigma_group"]][[n]]
      )
    },
    args = list(
      niter = niter,
      dat = dat,
      groups = igroups,
      mu0_global = mu0_global,
      Sigma0_global = Sigma0_global,
      mu0_group = mu0_group,
      Sigma0_group_inv = Sigma0_group_inv,
      v0_global = v0_global,
      S0_global = S0_global,
      v0_group = v0_group,
      S0_group = S0_group,
      setup_bygroup = setup_bygroup,
      progress = progress
    )
  )

  message("Running sampler...")
  raw_samples <- run_until_converged(
    sampler = sampler,
    model_type = "hier",
    inits = inits,
    nchains = nchains,
    max_attempts = max_attempts,
    save_progress = save_progress,
    threshold = threshold,
    keep_samples = keep_samples,
    autofit = autofit
  )

  message("Calculating correlation matrices...")
  raw_samples_corr <- add_correlations(raw_samples,
                                       hier = TRUE, ngroups = ngroup)

  message("Converting samples to coda mcmc.list object...")
  samples_mcmc <- results2mcmclist(raw_samples_corr, type = "hier")

  niter <- coda::niter(samples_mcmc)

  message("Preparing summary table...")
  summary_table <- summary_df(
    window(samples_mcmc, start = floor(niter / 2)),
    group = TRUE
  )

  stats <- c("Mean", "2.5%", "97.5%")
  mu_global_stats <- sapply(
    stats,
    function(x) summary2vec(summary_table, x,
                            variable == "mu", group == "global"),
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  Sigma_global_stats <- sapply(
    stats,
    function(x) summary2mat(summary_table, x,
                            variable == "Sigma", group == "global"),
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  Corr_global_stats <- sapply(
    stats,
    function(x) summary2mat(summary_table, x,
                            variable == "Corr", group == "global"),
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  get_mu_group <- function(grp) {
    sapply(
      stats,
      function(x) summary2vec(summary_table, x, variable == "mu", group == grp),
      simplify = FALSE,
      USE.NAMES = TRUE
    )
  }

  get_mat_group <- function(grp, var) {
    sapply(
      stats,
      function(x) summary2mat(summary_table, x, variable == var, group == grp),
      simplify = FALSE,
      USE.NAMES = TRUE
    )
  }

  mu_group_stats <- sapply(
    group_names,
    get_mu_group,
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  Sigma_group_stats <- sapply(
    group_names,
    get_mat_group,
    var = "Sigma",
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  Corr_group_stats <- sapply(
    group_names,
    get_mat_group,
    var = "Corr",
    simplify = FALSE,
    USE.NAMES = TRUE
  )

  list(
    summary_table = summary_table,
    stats = list(
      mu_global = mu_global_stats,
      Sigma_global = Sigma_global_stats,
      Corr_global = Corr_global_stats,
      mu_group = mu_group_stats,
      Sigma_group = Sigma_group_stats,
      Corr_group = Corr_group_stats
    ),
    samples = samples_mcmc
  )
}
