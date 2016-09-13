#' @import rstan
#' @export
customSTAN <- function(model_code, 
                       input_data, 
                       pars,
                       chains = 3,
                       model_name = "testmodel",
                       max.attempts = 50,
                       n_target = 5000,
                       rhat_max = 1.05,
                       iter = ceiling(n_target / chains),
                       ...){

    # STAN options for parallelism
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
    message("Running model...")
    continue <- TRUE
    attempt <- 1
    stanmodel <- stan_model(model_code = model_code,
                            model_name = model_name)
    while (continue & attempt <= max.attempts) {
        result <- sampling(stanmodel,
                           data = input_data,
                           iter = iter,
                           chains = chains,
                           pars = pars)
        stopifnot(result@mode == 0)
        result_summary <- summary(result)$summary
        # Ignore correlation coefficients
        # If covariance terms converge, correlations necessarily should also
        result_summary <- result_summary[!(grepl("Omega", rownames(result_summary))),]
        rhat <- result_summary[,"Rhat"]
        neff <- result_summary[,"n_eff"]
        lrhat <- rhat < rhat_max
        lneff <- neff >= n_target
        if (all(lrhat) & all(lneff)){
            message("Converged!")
            continue <- FALSE
        } else {
            if (!all(lrhat)) {
                message("The following values have not converged:")
                print(rhat[!lrhat])
            } else if (!all(lneff)) {
                message("The following values have insufficient samples")
                print(neff[!lneff])
            }
            continue <- TRUE
            iter <- round(iter * 2)
            message("Starting over with ", iter, " iterations")
        }
    } 
    if (attempt > max.attempts) {
        warning("Convergence was NOT achieved. Returning last result.")
    }
    return(result)
}

