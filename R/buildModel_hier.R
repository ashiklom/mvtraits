buildModel_hier <- function(dat, custom_inputs = list()) {

    pftcol <- grep("pft", colnames(dat))
    dat_all_missing <- apply(dat[,-pftcol], 1, function(x) all(is.na(x)))
    if (any(dat_all_missing)){
        warning("Some rows were missing entirely. Omitting these rows from analysis")
        dat <- dat[!dat_all_missing,]
    }
    datnp <- dat[,-pftcol]
    pftvec <- dat[,pftcol]

    dat_pres <- apply(datnp, 1, function(x) which(!is.na(x)))
    dat_pres_pft <- lapply(seq_along(dat_pres), 
                           function(x) list(ind = dat_pres[[x]],
                                            pft = pftvec[x]))
    namestring <- function(x) {
        sprintf("%s_pft_%02d",
                paste0(x[["ind"]], collapse = ""),
                x[["pft"]])
    }
    dat_pres_str <- sapply(dat_pres_pft, namestring)
    pres <- unique(dat_pres_str)
    pres_tab <- table(dat_pres_str)

    dat_pres_uniq_all <- unique(dat_pres_pft)
    dat_pres_uniq <- lapply(dat_pres_uniq_all, "[[", "ind")
    names(dat_pres_uniq) <- pres
    Npres <- sapply(dat_pres_uniq, length)
    pres <- pres[Npres > 0]

    pft_inds <- sapply(dat_pres_uniq_all, "[[", "pft")

    uni_ind <- which(Npres == 1)
    multi_ind <- which(Npres > 1)
    pres_uni <- pres[uni_ind]
    pres_multi <- pres[multi_ind]
    Npres_uni <- Npres[uni_ind]
    Npres_multi <- Npres[multi_ind]
    pft_inds_uni <- pft_inds[uni_ind]
    pft_inds_multi <- pft_inds[multi_ind]

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

    mu_definitions_uni <- sprintf("%s = mu_pft[%d,%d];", 
                                  mu_names_uni, 
                                  pft_inds_uni,
                                  ind_values_uni)
    mu_definitions_multi <- sprintf("%s = mu_pft[%d,%s];",
                                    mu_names_multi, 
                                    pft_inds_multi,
                                    ind_names_multi)
    mu_definitions <- c(mu_definitions_uni, mu_definitions_multi)

    Sigma_names_uni <- sprintf("sigma_%s", pres_uni)
    Sigma_names_multi <- sprintf("L_Sigma_%s", pres_multi)
    Sigma_declarations_uni <- sprintf("real %s;", Sigma_names_uni)
    Sigma_declarations_multi <- sprintf("matrix[%1$s,%1$s] %2$s;", 
                                        Npres_multi, Sigma_names_multi)

    Sigma_definitions_uni <- sprintf("%1$s = Sigma_pft[%3$d,%2$d,%2$d];",
                                     Sigma_names_uni, 
                                     ind_values_uni,
                                     pft_inds_uni)
    Sigma_definitions_multi <- sprintf("%1$s = cholesky_decompose(Sigma_pft[%3$d,%2$s,%2$s]);",
                                       Sigma_names_multi,
                                       ind_names_multi,
                                       pft_inds_multi)
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
    int<lower=0> Npft;

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
    vector[Npar] mu_global;
    cholesky_factor_corr[Npar] L_Omega_global;
    vector<lower=0>[Npar] sigma_vec_global;

    vector[Npar] mu_pft[Npft];
    cholesky_factor_corr[Npar] L_Omega_pft[Npft];
    vector<lower=0>[Npar] sigma_vec_pft[Npft];
}

transformed parameters {
    corr_matrix[Npar] Omega_global;
    cov_matrix[Npar] Sigma_global;

    corr_matrix[Npar] Omega_pft[Npft];
    cov_matrix[Npar] Sigma_pft[Npft];

    Omega_global = multiply_lower_tri_self_transpose(L_Omega_global);
    Sigma_global = quad_form_diag(Omega_global, sigma_vec_global);

    for (i in 1:Npft) {
        Omega_pft[i] = multiply_lower_tri_self_transpose(L_Omega_pft[i]);
        Sigma_pft[i] = quad_form_diag(Omega_pft[i], sigma_vec_pft[i]);
    }
}

model {
    matrix [Npar, Npar] L_Sigma_global;
    matrix [Npar, Npar] L_Sigma_pft[Npft];
    ",
    mu_declarations_uni,
    mu_declarations_multi,
    Sigma_declarations_uni,
    Sigma_declarations_multi,
    "

    // Prior
    mu_global ~ multi_normal(mu0, Sigma0);

    L_Omega_global ~ lkj_corr_cholesky(lkj_eta);
    sigma_vec_global ~ cauchy(cauchy_location, cauchy_scale);
    L_Sigma_global = diag_pre_multiply(sigma_vec_global, L_Omega_global);

    mu_pft ~ multi_normal_cholesky(mu_global, L_Sigma_global);

    for (i in 1:Npft) {
        L_Omega_pft[i] ~ lkj_corr_cholesky(lkj_eta);
        sigma_vec_pft[i] ~ cauchy(cauchy_location, cauchy_scale);
        L_Sigma_pft[i] = diag_pre_multiply(sigma_vec_pft[i], L_Omega_pft[i]);
    }
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

    Npar <- ncol(datnp)
    Npft <- length(unique(pftvec))
    default_inputs <- list(Npar = Npar,
                           Npft = Npft,
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
