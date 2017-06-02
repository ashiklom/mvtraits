#' @export
runModel <- function(model_type, 
                     dat,
                     groups = NA,
                     model_name = paste0('testmodel_', model_type),
                     chains = 3,
                     custom_inputs = list(),
                     ...){

    if (!model_type %in% c("uni", "multi", "hier")) {
        stop("Invalid model type. Must be 'uni', 'multi', or 'hier'")
    }

    if (model_type == 'hier' && is.na(groups)) {
        stop('If running hierarchical analysis, must supply groups vector.')
    }

    modbuild <- buildModel(model_type = model_type, 
                           dat = dat, 
                           groups = groups, 
                           custom_inputs = custom_inputs)
    model_code <- modbuild$model_code
    model_data <- modbuild$model_data

    variable_names <- switch(model_type, 
                             uni = c("mu", "sigma2"),
                             multi = c("mu", "Sigma", "Omega"),
                             hier = c("mu_group", "mu_global", 
                                      "Sigma_group", "Sigma_global", 
                                      "Omega_group", "Omega_global"))

    # Run model
    print(Sys.time())
    out <- try(customSTAN(model_code = model_code, 
                          model_data = model_data, 
                          pars = variable_names,
                          chains = chains,
                          model_name = model_name,
                          ...))
    return(out)
}

