#' @import data.table
#' @export
runModel <- function(model_type, 
                     try_data, 
                     pft_number, 
                     chains = 3,
                     ...){

    if (!model_type %in% c("uni", "multi", "hier")) {
        stop("Invalid model type. Must be 'uni', 'multi', or 'hier'")
    }

    if (!is.data.table(try_data)) try_data <- as.data.table(try_data)
    try_data[, pft := as.numeric(pft)]
    if (!is.na(pft_number)) {
        try_data <- try_data[pft == pft_number]
    } else {
        pft_number <- 0
    }
    #dat <- try_data[, !"pft", with=FALSE] %>% as.matrix()
    dat <- as.matrix(try_data, 
                     nrow = nrow(try_data),
                     ncol=ncol(try_data))
    model_name <- sprintf("%s_%02d", model_type, pft_number)

    #pftvec <- try_data[,pft]
    modbuild <- buildModel(model_type, dat, model_name)
    model_code <- modbuild$model_code
    model_data <- modbuild$model_data

    variable_names <- switch(model_type, 
                             uni = c("mu", "sigma2"),
                             multi = c("mu", "Sigma", "Omega"),
                             hier = c("mu_pft", "mu_global", 
                                      "Sigma_pft", "Sigma_global", 
                                      "Omega_pft", "Omega_global"))

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

