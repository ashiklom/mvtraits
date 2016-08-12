runModel <- function(model_type, try_data, pft_number,  n.chains = 3,
                     tau_obvs_miss = 0.01, tau_obvs_pres = 1000){
    if (!model_type %in% c("uni", "multi", "hier")) {
        stop("Invalid model type. Must be 'uni', 'multi', or 'hier'")
    }

    # Calculate global trait means for initial conditions
    global_means <- as.numeric(try_data[, lapply(.SD, mean, na.rm=TRUE),
                               .SDcols = traits])
    names(global_means) <- traits

    # Defining Gamma & Wishart parameters here so we can experiment
    # See test.wishart.R for more.
    # dgamma(gamma.shape,gamma.rate)
    # dwish(Wishart.rate,Wishart.df)

    n_traits <- 5
    n <- 1
    Wishart.rate <- diag(n, n_traits)
    Wishart.df <- n_traits
    mean <- n * Wishart.df
    gamma.shape <- Wishart.df/2 
    gamma.rate <- n/2   

    # Set up output storage directory
    out.dir <- paste0("output.n",n)
    if(!dir.exists(out.dir)) dir.create(out.dir)

    try_data[, pft := as.numeric(pft)]
    if (!is.na(pft_number)) {
        try_data <- try_data[pft == pft_number]
    }


    model <- system.file("models", paste0(model_type, ".bug"),
                         package = "mvtraits")
    print("Model path:")
    print(model)

    obvs <- try_data[, traits, with=FALSE] %>% as.matrix()
    n_obvs <- nrow(obvs)
    n_traits <- length(traits)

    data <- list(obvs = obvs,
                 n_traits = n_traits,
                 n_obvs = n_obvs)

    # Calculate trait mean values for iniital conditions
    trait_means <- try_data[, lapply(.SD, mean, na.rm = T), .SDcols = traits
                        ][, lapply(.SD, nan2na)] %>%
        as.numeric() %>%
        replace.with.global(global_means)

    if (model_type == "uni") {
        # Univariate
        data <- modifyList(data, list(gamma.shape = gamma.shape,
                                      gamma.rate = gamma.rate))       # scale/rate = n
        inits = function() list(mu_trait = as.numeric(trait_means))
        variable_names <- c("mu_trait", "sigma2_obvs")

    } else if (model_type %in% c("multi", "hier")) {
        # Missing values
        miss <- which(is.na(obvs), arr.ind = TRUE)
        pres <- which(!is.na(obvs), arr.ind = TRUE)
        n_miss <- nrow(miss)
        n_pres <- nrow(pres)
        data <- modifyList(data, 
                           list(miss = miss,
                                pres = pres,
                                n_miss = n_miss,
                                n_pres = n_pres,
                                mu0 = rep(0,n_traits), 
                                Sigma0 = diag(0.001,n_traits),
                                Wishart.rate = Wishart.rate,
                                Wishart.df = Wishart.df,
                                tau_obvs_miss = tau_obvs_miss,
                                tau_obvs_pres = tau_obvs_pres))

        if (model_type == "multi") {
            # Multivariate
            inits = function() list(mu_trait = as.numeric(trait_means))
            variable_names <- c("mu_trait", "Sigma_trait")

        } else if (model_type == "hier") {
            # Hierarchical
            pft_means <- try_data[, lapply(.SD, mean, na.rm = TRUE), 
                            by = pft, .SDcols = "traits"
                            ][, lapply(.SD, nan2na)
                            ][, traits, with = FALSE] %>%
                        as.matrix() %>% 
                        apply(1, replace.with.global, global_means) %>%
                        t()
            n_pfts <- nrow(pft_means)
            data <- modifyList(data, list(pft_obvs = try_data[,pft]))
            inits = function() list(mu_trait = as.numeric(trait_means),
                                    mu_pft_trait = as.matrix(pft_means))
            variable_names <- c("mu_trait","Sigma_trait", 
                                "mu_pft_trait", "Sigma_pft")
        }
    }

    # Run model
    print(Sys.time())
    out <- try(customJAGS(model = model, data = data, inits = inits,
                          n.chains = n.chains,
                          variable.names = variable_names))
    return(out)
}

