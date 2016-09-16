buildModel_uni <- function(dat, custom_inputs = list()) {

    if ("pft" %in% colnames(dat)) {
        dat <- dat[, colnames(dat) != "pft", drop=FALSE]
    }

    dat_all_missing <- apply(dat, 1, function(x) all(is.na(x)))
    if (any(dat_all_missing)){
        warning("Some rows were missing entirely. Omitting these rows from analysis")
        dat <- dat[!dat_all_missing,]
    }

    # Separate each column into a list item
    dat_cols_full <- split(dat, rep(1:ncol(dat), each = nrow(dat)))

    # Remove NA values from each column
    colind <- 1:ncol(dat)
    dat_cols <- lapply(dat_cols_full, function(x) x[!is.na(x)])
    names(dat_cols) <- sprintf("dat_%d", colind)

    # Indices of columns that have at least one observation
    dat_lengths_all <- sapply(dat_cols, length)
    use <- dat_lengths_all > 0
    dat_lengths <- dat_lengths_all[use]

    coluse <- colind[use]

    # Coerce every column to an array to keep STAN happy
    dat_in <- mapply(array, dat_cols[use], dat_lengths)

    dat_names <- names(dat_in)
    dat_declarations <- sprintf("real %s[%d];", dat_names, dat_lengths)

    mu_names <- sprintf("mu_%d", coluse)
    mu_declarations <- sprintf("real %s;", mu_names)
    mu_definitions <- sprintf("%s = mu[%d];", mu_names, coluse)
    
    sigma_names <- sprintf("sigma_%d", coluse)
    sigma_declarations <- sprintf("real %s;", sigma_names)
    sigma_definitions <- sprintf("%s = sigma[%d];", sigma_names, coluse)

    sampling_statements <- sprintf("%s ~ normal(%s, %s);",
                                   dat_names, mu_names, sigma_names)
    
    full_model_string <- c(
" data {
    int<lower=0> Npar;

    vector[Npar] mu0;
    vector[Npar] sigma0;

    real<lower=0> cauchy_location;
    real<lower=0> cauchy_scale;
",

    dat_declarations,
"}

parameters {
    vector[Npar] mu;
    vector<lower=0>[Npar] sigma;
}

model {
    ",
    mu_declarations,
    sigma_declarations,
    "

    // Prior
    mu ~ normal(mu0, sigma0);
    sigma ~ cauchy(cauchy_location, cauchy_scale);
    ",

    mu_definitions,
    sigma_definitions,
    sampling_statements,
    "
}

generated quantities {
    vector[Npar] sigma2;
    sigma2 = sigma .* sigma;
}
")

    model_code <- paste(full_model_string, collapse = "\n")

    Npar <- ncol(dat)
    default_inputs <- list(Npar = Npar,
                        mu0 = rep(0, Npar),
                        sigma0 = rep(sqrt(1000), Npar),
                        cauchy_location = 0,
                        cauchy_scale = 2.5)

    common_inputs <- modifyList(default_inputs, custom_inputs)

    model_data <- c(common_inputs, dat_in)

    out <- list(model_code = model_code, model_data = model_data)
    return(out)
}
