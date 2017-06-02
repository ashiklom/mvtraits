buildModel_hier <- function(dat, groups, custom_inputs = list()) {

    Ngroup <- dplyr::n_distinct(groups)
    params <- colnames(dat)
    if (is.null(params)) {
        params <- paste0('V', 1:ncol(dat))
    }

    stopifnot(length(groups) == nrow(dat), Ngroup == max(groups))

    dat_df <- dplyr::as_data_frame(dat)

    dat_df_full <- dplyr::mutate(dat_df, 
                                 group = groups,
                                 rn = 1:nrow(dat_df))

    dat_long <- tidyr::gather(dat_df_full, variable, value, -rn, -group)

    dat_nested <- tidyr::nest(dat_long, variable, value)

    dat_present <- dplyr::mutate(dat_nested, present = purrr::map(data, ~which(is.na(.$value))),
                                 present_str = purrr::map_chr(present, paste, collapse = "_"),
                                 npresent = purrr::map_int(present, length),
                                 present_str_group = sprintf('g_%d_p_%s', group, present_str)
                                 )

    dat_present_sub <- dplyr::filter(dplyr::select(dat_present, -data), npresent > 0)
    diff_rows <- nrow(dat_present) - nrow(dat_present_sub)
    if (diff_rows > 0) {
        warning(diff_rows, ' rows dropped because all values were missing.')
    }

    dat_merge <- dplyr::left_join(dat_present_sub, dat_df_full)
    #dat_renested <- tidyr::nest(dat_merge, -rn, -present, -group, -present_str_group, -npresent)
    dat_renested <- tidyr::nest(dplyr::select(dat_merge, -rn, -present), dplyr::one_of(params))

    dat_processed <- 
        dplyr::mutate(dat_renested,
                      present = purrr::map(present_str, function(x) as.integer(unlist(strsplit(x, "_")))),
                      Nrow_str = sprintf("Nrow_%s", present_str_group),
                      Nrow_dec = sprintf("int<lower=0> %s;", Nrow_str),
                      Nrow_values = purrr::map_dbl(data, nrow),
                      dat_str = sprintf("dat_%s", present_str_group),
                      dat_dec = dplyr::if_else(npresent == 1, 
                                               sprintf('real %s[%s];', dat_str, Nrow_str),
                                               sprintf('vector[%d] %s[%s];', npresent, dat_str, Nrow_str)),
                      ind_str = sprintf('ind_%s', present_str_group),
                      ind_dec = sprintf('int %s[%s];', ind_str, Nrow_str),
                      mu_str = sprintf('mu_%s', present_str_group),
                      mu_dec = dplyr::if_else(npresent == 1,
                                              sprintf('real %s;', mu_str),
                                              sprintf('vector[%d] %s;', npresent, mu_str)),
                      mu_def = sprintf('%s = mu_group[%d,%s];', mu_str, group, ind_str),
                      Sigma_str = dplyr::if_else(npresent == 1,
                                                 sprintf('sigma_vec_%s', present_str_group),
                                                 sprintf('L_Sigma_%s', present_str_group)),
                      Sigma_dec = dplyr::if_else(npresent == 1,
                                                 sprintf('real %s;', Sigma_str),
                                                 sprintf('matrix[%1$s,%1$s] %2$s;', npresent, Sigma_str)),
                      Sigma_def = dplyr::if_else(npresent == 1,
                                                 sprintf('%1$s = Sigma_group[%3$d,%2$s,%2$s];',
                                                         Sigma_str, ind_str, group),
                                                 sprintf('%1$s = cholesky_decompose(Sigma_group[%3$d,%2$s,%2$s,%2$s]);',
                                                         Sigma_str, ind_str, group)),
                      sample_statement = dplyr::if_else(npresent == 1,
                                                        sprintf('%s ~ normal(%s, %s);', 
                                                                dat_str, mu_str, Sigma_str),
                                                        sprintf('%s ~ multi_normal_cholesky(%s, %s);',
                                                                dat_str, mu_str, Sigma_str))
                      )
    #glimpse(dat_processed)

    model_dat <- dplyr::arrange(dat_processed, npresent, group, present_str_group)

    full_model_string <- c(
                           " data {
                           int<lower=0> Npar;
                           int<lower=0> Ngroup;

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
                               vector[Npar] mu_global;
                               cholesky_factor_corr[Npar] L_Omega_global;
                               vector<lower=0>[Npar] sigma_vec_global;

                               vector[Npar] mu_group[Ngroup];
                               cholesky_factor_corr[Npar] L_Omega_group[Ngroup];
                               vector<lower=0>[Npar] sigma_vec_group[Ngroup];
                           }

                           transformed parameters {
                               corr_matrix[Npar] Omega_global;
                               cov_matrix[Npar] Sigma_global;

                               corr_matrix[Npar] Omega_group[Ngroup];
                               cov_matrix[Npar] Sigma_group[Ngroup];

                               Omega_global = multiply_lower_tri_self_transpose(L_Omega_global);
                               Sigma_global = quad_form_diag(Omega_global, sigma_vec_global);

                               for (i in 1:Ngroup) {
                                   Omega_group[i] = multiply_lower_tri_self_transpose(L_Omega_group[i]);
                                   Sigma_group[i] = quad_form_diag(Omega_group[i], sigma_vec_group[i]);
                               }
                           }

                           model {
                               matrix [Npar, Npar] L_Sigma_global;
                               matrix [Npar, Npar] L_Sigma_group[Ngroup];
                               ",
                               unique(model_dat[['mu_dec']]),
                               unique(model_dat[['Sigma_dec']]),
                               "

                               // Prior
                               mu_global ~ multi_normal(mu0, Sigma0);

                               L_Omega_global ~ lkj_corr_cholesky(lkj_eta);
                               sigma_vec_global ~ cauchy(cauchy_location, cauchy_scale);
                               L_Sigma_global = diag_pre_multiply(sigma_vec_global, L_Omega_global);

                               mu_group ~ multi_normal_cholesky(mu_global, L_Sigma_global);

                               for (i in 1:Ngroup) {
                                   L_Omega_group[i] ~ lkj_corr_cholesky(lkj_eta);
                                   sigma_vec_group[i] ~ cauchy(cauchy_location, cauchy_scale);
                                   L_Sigma_group[i] = diag_pre_multiply(sigma_vec_group[i], L_Omega_group[i]);
                               }
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
                                                  Ngroup = Ngroup,
                                                  mu0 = rep(0, Npar),
                                                  Sigma0 = diag(1000, Npar),
                                                  cauchy_location = 0,
                                                  cauchy_scale = 2.5,
                                                  lkj_eta = 1)
                           common_inputs <- modifyList(default_inputs, custom_inputs)

                           Nrow_in <- model_dat[['npresent']]
                           names(Nrow_in) <- model_dat[['Nrow_str']]

                           dat_in <- model_dat[['data']]
                           names(dat_in) <- model_dat[['dat_str']]

                           ind_in <- model_dat[['present']]
                           names(ind_in) <- model_dat[['ind_str']]

                           return_model_data <- c(common_inputs, Nrow_in, dat_in, ind_in)
                           out <- list(model_code = model_code, model_data = return_model_data)

                           return(out)

}
