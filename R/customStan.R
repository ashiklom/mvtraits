#' @import rstan
#' @export
customSTAN <- function(model_code, 
                       model_data, 
                       pars,
                       chains = 3,
                       model_name = "testmodel",
                       max.attempts = 10,
                       n_target = 1,
                       rhat_max = 1.05,
                       iter = 5000,
                       save_each = TRUE,
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
                           data = model_data,
                           iter = iter,
                           chains = chains,
                           pars = pars,
                           ...)
        if (result@mode != 0) {
            #print(result@stanmodel)
            stop("Error in stan model")
        }
        result_summary <- summary(result)$summary
        if (save_each) {
            now_time <- strftime(Sys.time(), "%Y_%m_%d_%H_%M")
            tempfilename <- sprintf("%s.%s.rds", model_name, now_time)
            saveRDS(result, file = tempfilename)
        }
        # Ignore correlation coefficients
        # If covariance terms converge, correlations must also
        result_summary <- result_summary[!(grepl("Omega|lp__", rownames(result_summary))),]
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

