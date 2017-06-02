buildModel_multi <- function(dat, custom_inputs = list()) {

    stopifnot(length(dim(dat)) > 1)

    params <- colnames(dat)
    if (is.null(params)) {
        params <- paste0('V', 1:ncol(dat))
    }

    dat_df <- dplyr::as_data_frame(dat)

    dat_df_full <- dplyr::mutate(dat_df, rn = 1:nrow(dat_df))

    dat_long <- tidyr::gather(dat_df_full, variable, value, -rn)

    dat_nested <- tidyr::nest(dat_long, variable, value)

    dat_present <- dplyr::mutate(dat_nested, 
                                 present = purrr::map(data, ~which(!is.na(.$value))),
                                 present_str = purrr::map_chr(present, paste, collapse = "_"),
                                 npresent = purrr::map_int(present, length)
                                 )

    dat_present_sub <- dplyr::filter(dplyr::select(dat_present, -data), npresent > 0)
    diff_rows <- nrow(dat_present) - nrow(dat_present_sub)
    if (diff_rows > 0) {
        warning(diff_rows, ' rows dropped because all values were missing.')
    }

    dat_merge <- dplyr::left_join(dat_present_sub, dat_df_full)
    dat_renested <- tidyr::nest(dplyr::select(dat_merge, -rn, -present), dplyr::one_of(params))

    matrify <- function(x, y) {
        ymat <- as.matrix(y)
        drp <- TRUE
        if (nrow(ymat) == 1) drp <- FALSE
        if (length(x) == 1) drp <- TRUE
        d <- ymat[, x, drop = drp]
        return(as.array(d))
    }

    dat_processed <- 
        dplyr::mutate(dat_renested,
                      present = purrr::map(present_str, function(x) as.array(as.integer(unlist(strsplit(x, "_"))))),
                      data_sub = purrr::map2(present, data, matrify),
                      Nrow_str = sprintf("Nrow_%s", present_str),
                      Nrow_dec = sprintf("int<lower=0> %s;", Nrow_str),
                      Nrow_values = purrr::map_dbl(data, nrow),
                      dat_str = sprintf("dat_%s", present_str),
                      dat_dec = dplyr::if_else(npresent == 1, 
                                               sprintf('real %s[%s];', dat_str, Nrow_str),
                                               sprintf('vector[%d] %s[%s];', npresent, dat_str, Nrow_str)),
                      ind_str = sprintf('ind_%s', present_str),
                      ind_dec = sprintf('int %s[%d];', ind_str, npresent),
                      mu_str = sprintf('mu_%s', present_str),
                      mu_dec = dplyr::if_else(npresent == 1,
                                              sprintf('real %s;', mu_str),
                                              sprintf('vector[%d] %s;', npresent, mu_str)),
                      mu_def = dplyr::if_else(npresent == 1,
                                              sprintf('%s = mu[%d];', mu_str,
                                                      purrr::map_int(present, 1)),
                                              sprintf('%s = mu[%s];', mu_str, ind_str)),
                      Sigma_str = dplyr::if_else(npresent == 1,
                                                 sprintf('sigma_vec_%s', present_str),
                                                 sprintf('L_Sigma_%s', present_str)),
                      Sigma_dec = dplyr::if_else(npresent == 1,
                                                 sprintf('real %s;', Sigma_str),
                                                 sprintf('matrix[%1$s,%1$s] %2$s;', npresent, Sigma_str)),
                      Sigma_def = dplyr::if_else(npresent == 1,
                                                 sprintf('%1$s = Sigma[%2$d,%2$d];',
                                                         Sigma_str, purrr::map_int(present, 1)),
                                                 sprintf('%1$s = cholesky_decompose(Sigma[%2$s,%2$s]);',
                                                         Sigma_str, ind_str)),
                      sample_statement = dplyr::if_else(npresent == 1,
                                                        sprintf('%s ~ normal(%s, %s);', 
                                                                dat_str, mu_str, Sigma_str),
                                                        sprintf('%s ~ multi_normal_cholesky(%s, %s);',
                                                                dat_str, mu_str, Sigma_str))
                      )
    #glimpse(dat_processed)

    model_dat <- dplyr::arrange(dat_processed, npresent, present_str)

    full_model_string <- c(
                           " data {
                           int<lower=0> Npar;

                           vector[Npar] mu0;
                           cov_matrix[Npar] Sigma0;

                           real<lower=0> cauchy_location;
                           real<lower=0> cauchy_scale;
                           real<lower=0> lkj_eta;
                           ",
                           unique(model_dat[['Nrow_dec']]),
                           unique(model_dat[['dat_dec']]),
                           unique(model_dat[['ind_dec']]),
                           "}

                           parameters {
                               vector[Npar] mu;
                               cholesky_factor_corr[Npar] L_Omega;
                               vector<lower=0>[Npar] sigma_vec;
                           }

                           transformed parameters {
                               corr_matrix[Npar] Omega;
                               cov_matrix[Npar] Sigma;

                               Omega = multiply_lower_tri_self_transpose(L_Omega);
                               Sigma = quad_form_diag(Omega, sigma_vec);
                           }

                           model {
                               matrix [Npar, Npar] L_Sigma;
                               ",
                               unique(model_dat[['mu_dec']]),
                               unique(model_dat[['Sigma_dec']]),
                               "

                               // Prior
                               mu ~ multi_normal(mu0, Sigma0);

                               L_Omega ~ lkj_corr_cholesky(lkj_eta);
                               sigma_vec ~ cauchy(cauchy_location, cauchy_scale);
                               L_Sigma = diag_pre_multiply(sigma_vec, L_Omega);

                               ",
                               unique(model_dat[['mu_def']]),
                               unique(model_dat[['Sigma_def']]),
                               unique(model_dat[['sample_statement']]),
                               "
                           }
                           ")

                           model_code <- paste(full_model_string, collapse = '\n')

                           Npar <- ncol(dat)
                           default_inputs <- list(Npar = Npar,
                                                  mu0 = rep(0, Npar),
                                                  Sigma0 = diag(1000, Npar),
                                                  cauchy_location = 0,
                                                  cauchy_scale = 2.5,
                                                  lkj_eta = 1)
                           common_inputs <- modifyList(default_inputs, custom_inputs)

                           Nrow_in <- model_dat[['Nrow_values']]
                           names(Nrow_in) <- model_dat[['Nrow_str']]

                           dat_in <- model_dat[['data_sub']]
                           names(dat_in) <- model_dat[['dat_str']]

                           ind_in <- model_dat[['present']]
                           names(ind_in) <- model_dat[['ind_str']]

                           return_model_data <- c(common_inputs, Nrow_in, dat_in, ind_in)
                           out <- list(model_code = model_code, model_data = return_model_data)

                           return(out)

}
