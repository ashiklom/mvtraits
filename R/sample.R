sample_mvnorm <- function(niter, dat,
                          mu, Sigma,
                          mu0, Sigma0_inv,
                          v0, S0,
                          setup,
                          progress = FALSE) {
  # Re-arrange dimensions of Sigma_samp to fit C code
  ## result <- r_sample_mvnorm(niter, dat, mu, Sigma,
  ##                           mu0, Sigma0_inv, v0, S0,
  ##                           setup, progress = progress)
  # Capture annoying "inv_sympd" warnings
  sampler_stderr <- capture.output({
    result <- c_sample_mvnorm(niter, dat, mu, Sigma,
                              mu0, Sigma0_inv, v0, S0,
                              setup, progress)
  }, type = "message")
  ## attr(result, "stderr") <- sampler_stderr
  # Add names to results
  params <- names(mu)
  colnames(result$mu) <- params
  colnames(result$Sigma) <- lowertri_names(params)
  return(result)
}

sample_mvnorm_hier <- function(niter, dat, groups,
                               mu_global, Sigma_global,
                               mu_group, Sigma_group,
                               mu0_global, Sigma0_global_inv,
                               mu0_group, Sigma0_group_inv,
                               v0_global, S0_global,
                               v0_group, S0_group,
                               setup_bygroup,
                               progress = FALSE) {

  # Convert dimensions
  c_Sigma_group <- aperm(Sigma_group, c(2, 3, 1))
  c_Sigma0_group_inv <- aperm(Sigma0_group_inv, c(2, 3, 1))
  c_S0_group <- aperm(S0_group, c(2, 3, 1))

  # Run sampler
  sampler_stderr <- capture.output({
    result <- c_sample_mvnorm_hier(niter, dat, groups,
                                   mu_global, Sigma_global,
                                   mu_group, c_Sigma_group,
                                   mu0_global, Sigma0_global_inv,
                                   mu0_group, c_Sigma0_group_inv,
                                   v0_global, S0_global,
                                   v0_group, c_S0_group,
                                   setup_bygroup,
                                   progress)
  }, type = "message")

  ## Add names
  params <- names(mu_global)
  nparams <- length(params)
  ugroups <- rownames(mu_group)
  ngroups <- length(ugroups)
  params_rep <- rep(params, each = ngroups)
  groups_rep <- rep(ugroups, nparams)
  params_groups <- paste(groups_rep, params_rep, sep = varsep)
  sigma_params <- lowertri_names(params, diag = TRUE)
  sigma_params_rep <- rep(sigma_params, ngroups)
  sigma_groups_rep <- rep(ugroups, each = length(sigma_params))
  sigma_params_groups <- paste(sigma_groups_rep, sigma_params_rep, sep = varsep)
  colnames(result$mu_global) <- params
  colnames(result$Sigma_global) <- sigma_params
  colnames(result$mu_group) <- params_groups
  colnames(result$Sigma_group) <- sigma_params_groups
  return(result)
}
