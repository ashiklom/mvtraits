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
    }
    dat <- try_data[, !"pft", with=FALSE] %>% as.matrix()

    # Define common input_data
    Nrow <- nrow(dat)
    Npar <- ncol(dat)
    pftvec <- try_data[,pft]
    Npft <- length(unique(pftvec))
    sigma_0 <- 1000
    input_data <- list(Nrow = Nrow,
                       Npar = Npar,
                       pftvec = pftvec,
                       Npft = Npft,
                       # Priors
                       mu0 = rep(0, Npar),
                       sigma0 = rep(sigma_0, Npar),     # uni
                       Sigma0 = diag(sigma_0, Npar),    # multi/hier
                       cauchy_location = 0,
                       cauchy_scale = 2.5,
                       lkj_eta = 1)

    model_file <- system.file("models", paste0(model_type, ".stan"),
                         package = "mvtraits")
    if (file.exists(model_file)) {
        message("Model path: ", model_file)
        model_code <- readChar(model_file, file.info(model_file)$size)
    } else {
        message("Model file not found. Assuming model_code string")
        model_code <- model
    }

    if (model_type %in% c("multi", "hier")) {
        # Ragged array storage
        # See STAN manual, section 13.2 (Ragged Data Structures), pg. 150-151 
        dat_present <- !is.na(dat)
        # Need the following:
        #   sizevec -- VECTOR of number of present traits per row [nrow]
        input_data$sizevec <- rowSums(dat_present)

        #   indvec -- VECTOR of positions of present traits [nrow x ncol - nmiss]
        indmat <- which(dat_present, arr.ind=TRUE)
        indmat <- indmat[order(indmat[,"row"], indmat[,"col"]),]
        input_data$indvec <- indmat[,"col"]

        #   datvec -- VECTOR of all of the values [nrow x ncol - nmiss]
        dat_t <- t(dat)
        datvec <- dat_t[!is.na(dat_t)]
        input_data$y <- datvec
        input_data$Ndata <- length(datvec)
    } else {
        input_data$pres <- which(!is.na(dat), arr.ind = TRUE)
        input_data$Ndata <- nrow(input_data$pres)
        dat[is.na(dat)] <- -9999   # Values ignored in model
        input_data$y <- dat
    }

    if (model_type == "uni") {
        # Univariate
        variable_names <- c("mu", "sigma2")
    } else if (model_type %in% c("multi", "hier")) {
        if (model_type == "multi") {
            # Multivariate
            variable_names <- c("mu", "Sigma", "Omega")

        } else if (model_type == "hier") {
            # Hierarchical
            variable_names <- c("mu_pft","Sigma_pft", "Omega_pft", 
                                "mu_global", "Sigma_global", "Omega_global")
        }
    }

    # Run model
    print(Sys.time())
    out <- try(customSTAN(model_code = model_code, 
                          input_data = input_data, 
                          pars = variable_names,
                          chains = chains,
                          ...))
    return(out)
}

