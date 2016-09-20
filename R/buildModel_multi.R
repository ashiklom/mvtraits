buildModel_multi <- function(dat, custom_inputs = list()) {

    if ("pft" %in% colnames(dat)) {
        dat <- dat[, colnames(dat) != "pft", drop=FALSE]
    }

    dat_all_missing <- apply(dat, 1, function(x) all(is.na(x)))
    if (any(dat_all_missing)){
        warning("Some rows were missing entirely. Omitting these rows from analysis")
        dat <- dat[!dat_all_missing,]
    }
    dat_pres <- plyr::alply(dat, 1, function(x) which(!is.na(x)))
    dat_pres_str <- sapply(dat_pres, paste0, collapse = "")
    pres <- unique(dat_pres_str)
    pres_tab <- table(dat_pres_str)

    dat_pres_uniq <- unique(dat_pres)
    names(dat_pres_uniq) <- pres
    Npres <- sapply(dat_pres_uniq, length)
    pres <- pres[Npres > 0]

    uni_ind <- which(Npres == 1)
    multi_ind <- which(Npres > 1)
    pres_uni <- pres[uni_ind]
    pres_multi <- pres[multi_ind]
    Npres_uni <- Npres[uni_ind]
    Npres_multi <- Npres[multi_ind]

    Nrow <- pres_tab
    Nrow_names <- sprintf("Nrow_%s", pres)
    Nrow_declarations <- sprintf("int<lower=0> %s;", Nrow_names)
    Nrow_in <- as.list(pres_tab)[pres]
    names(Nrow_in) <- Nrow_names

    dat_names_uni <- sprintf("dat_%s", pres_uni)
    dat_names_multi <- sprintf("dat_%s", pres_multi)
    # TODO: Maybe fix
    dat_declarations_uni <- sprintf("real %s[%s];", 
                                    dat_names_uni, Nrow_names[uni_ind])
    dat_declarations_multi <- sprintf("vector[%s] %s[%s];", 
                                      Npres_multi, dat_names_multi,
                                      Nrow_names[multi_ind])

    ind_values_uni <- unlist(dat_pres_uniq[uni_ind])
    ind_names_multi <- sprintf("ind_%s", pres_multi)
    ind_declarations <- sprintf("int %s[%d];", ind_names_multi, Npres_multi)


    mu_names_uni <- sprintf("mu_%s", pres_uni)
    mu_names_multi <- sprintf("mu_%s", pres_multi)
    mu_declarations_uni <- sprintf("real %s;", mu_names_uni)
    mu_declarations_multi <- sprintf("vector[%d] %s;",
                                     Npres_multi, mu_names_multi)

    mu_definitions_uni <- sprintf("%s = mu[%d];", mu_names_uni, ind_values_uni)
    mu_definitions_multi <- sprintf("%s = mu[%s];",
                                    mu_names_multi, ind_names_multi)
    mu_definitions <- c(mu_definitions_uni, mu_definitions_multi)

    Sigma_names_uni <- sprintf("sigma_%s", pres_uni)
    Sigma_names_multi <- sprintf("L_Sigma_%s", pres_multi)
    Sigma_declarations_uni <- sprintf("real %s;", Sigma_names_uni)
    Sigma_declarations_multi <- sprintf("matrix[%1$s, %1$s] %2$s;", 
                                        Npres_multi, Sigma_names_multi)

    Sigma_definitions_uni <- sprintf("%1$s = Sigma[%2$d,%2$d];",
                                     Sigma_names_uni, ind_values_uni)
    Sigma_definitions_multi <- sprintf("%1$s = cholesky_decompose(Sigma[%2$s,%2$s]);",
                                       Sigma_names_multi, ind_names_multi)
    Sigma_definitions <- c(Sigma_definitions_uni, Sigma_definitions_multi)

    sampling_statements_uni <- sprintf("%s ~ normal(%s, %s);",
                                       dat_names_uni, mu_names_uni,
                                       Sigma_names_uni)
    sampling_statements_multi <- 
        sprintf("%s ~ multi_normal_cholesky(%s, %s);",
                dat_names_multi, mu_names_multi, Sigma_names_multi)

    full_model_string <- c(
" data {
    int<lower=0> Npar;

    vector[Npar] mu0;
    cov_matrix[Npar] Sigma0;

    real<lower=0> cauchy_location;
    real<lower=0> cauchy_scale;
    real<lower=0> lkj_eta;
",

    Nrow_declarations,
    dat_declarations_uni,
    dat_declarations_multi,
    ind_declarations,
"}

parameters {
    vector[Npar] mu;
    vector<lower=0>[Npar] sigma_vec;
    cholesky_factor_corr[Npar] L_Omega;
}

transformed parameters {
    corr_matrix[Npar] Omega;
    cov_matrix[Npar] Sigma;

    Omega = multiply_lower_tri_self_transpose(L_Omega);
    Sigma = quad_form_diag(Omega, sigma_vec);
}

model {
",
    mu_declarations_uni,
    mu_declarations_multi,
    Sigma_declarations_uni,
    Sigma_declarations_multi,

"
    // Prior
    mu ~ multi_normal(mu0, Sigma0);

    L_Omega ~ lkj_corr_cholesky(lkj_eta);
    sigma_vec ~ cauchy(cauchy_location, cauchy_scale);

    ",

    mu_definitions,
    Sigma_definitions,
    sampling_statements_uni,
    sampling_statements_multi,
    "
}
")

    model_code <- paste(full_model_string, collapse = "\n")
    #cat(model_code)

    Npar <- ncol(dat)
    default_inputs <- list(Npar = Npar,
                        mu0 = rep(0, Npar),
                        Sigma0 = diag(1000, Npar),
                        cauchy_location = 0,
                        cauchy_scale = 2.5,
                        lkj_eta = 1)
    common_inputs <- modifyList(default_inputs, custom_inputs)

    select_dat <- function(x) {
        drp <- TRUE
        irow <- which(dat_pres_str == x)
        icol <- dat_pres_uniq[[x]]
        if (length(irow) == 1) drp <- FALSE
        if (length(icol) == 1) drp <- TRUE
        d <- dat[irow, icol, drop = drp]
        if (length(d) == 1) d <- array(d, 1)
        return(d)
    }

    dat_in <- lapply(pres, select_dat)
    names(dat_in) <- sprintf("dat_%s", pres)

    ind_in <- dat_pres_uniq[pres_multi]
    names(ind_in) <- ind_names_multi

    model_data <- c(common_inputs, Nrow_in, dat_in, ind_in)

    out <- list(model_code = model_code, model_data = model_data)
    return(out)
}
