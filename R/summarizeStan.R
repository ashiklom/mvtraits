#' @export
summarizeStan <- function(filename) {
    library(rstan)
    library(data.table)
    library(dplyr)
    library(tidyr)

    # Parse model string
    model_string <- ".*[/](uni|multi|hier)_?([[:digit:]]{0,2})\\.Rdata"
    model_type <- gsub(model_string, "\\1", filename)
    model_pft <- as.numeric(gsub(model_string, "\\2", filename))
    if (is.na(model_pft)) model_pft <- 0

    # Paramters list
    npar <- 5
    vpar <- 1:npar
    npft <- 35
    vpft <- 1:npft
    vpar_pft <- rep(vpar, npft)
    vpft_par <- rep(vpft, each = npar)
    mm <- matrix(numeric(), npar, npar)
    spar <- which(lower.tri(mm, diag=TRUE), arr.ind=TRUE)
    opar <- which(lower.tri(mm, diag=FALSE), arr.ind=TRUE)
    spar_pft <- rep(vpft, each = nrow(spar))
    opar_pft <- rep(vpft, each = nrow(opar))
    params <- list(mu = sprintf("mu[%d]", vpar),
                   mu_global = sprintf("mu_global[%d]", vpar),
                   mu_pft = sprintf("mu_global[%d,%d]", 
                                      vpft_par, vpar_pft),
                   sigma2 = sprintf("sigma2[%d]", vpar),
                   Sigma = sprintf("Sigma[%d,%d]", spar[,1], spar[,2]),
                   Omega = sprintf("Omega[%d,%d]", opar[,1], opar[,2]),
                   Sigma_global = sprintf("Sigma_global[%d,%d]",
                                          spar[,1], spar[,2]),
                   Omega_global = sprintf("Omega_global[%d,%d]",
                                          opar[,1], opar[,2]),
                   Sigma_pft = sprintf("Sigma_pft[%d,%d,%d]",
                                         spar_pft, spar[,1], spar[,2]),
                   Omega_pft = sprintf("Omega_pft[%d,%d,%d]",
                                         opar_pft, opar[,1], opar[,2]))

    # Load model and grab the relevant parameters
    load(filename)
    allpars <- out@model_pars
    outsum <- summary(out)$summary
    result <- as.data.table(outsum, keep.rownames = TRUE)
    result <- result[rn %in% unlist(params)]
    rm(out); gc()
    result[, c("model_type", "model_pft_num") := list(model_type, model_pft)]

    # Parse parameter type
    rnstring <- "^(.*)\\[([[:digit:]]{,2}),?([[:digit:]])?,?([[:digit:]])?\\]"
    result[, var_type := gsub(rnstring, "\\1", rn)]

    getnum <- function(num, string) {
        num <- num + 1
        as.numeric(gsub(rnstring, paste0("\\", num), string))
    }

    cols <- c("trait", "var_pft_num")
    sg <- c("Sigma", "Sigma_global")
    og <- c("Omega", "Omega_global")
    sog <- c(sg, og)
    sop <- c("Sigma_pft", "Omega_pft")

    result[,n1 := getnum(1,rn)][,n2 := getnum(2,rn)][,n3 := getnum(3,rn)]
    result[var_type %in% c("mu", "mu_global"), 
           c(cols) := list(traits_nolog[n1], 0)]

    ptrait <- function(n1, n2) paste(traits_nolog[n1], traits_nolog[n2], sep=".")

    # Assign Sigma and Omega trait names
    result[var_type %in% sog, c(cols) := list(ptrait(n1, n2), 0)]

    # Hierarchical indices
    result[var_type == "mu_pft", c(cols) := list(traits_nolog[n2], n1)]
    result <- result[!(var_type=="Sigma_pft" & n2 > n3)]
    result <- result[!(var_type=="Omega_pft" & n2 >= n3)]
    result[var_type %in% sop, c(cols) := list(ptrait(n2, n3), n1)]
    
    result <- result[!is.na(trait)]

    # Assign PFT name
    result[model_type %in% c("uni", "multi"), pft_num := model_pft_num]
    result[model_type == "hier", pft_num := var_pft_num]
    result[pft_num == 0, PFT := "global"]
    result[pft_num != 0, PFT := pft.names[pft_num]]

    # Fix column names and remove unneeded columns
    names_dict <- c("model_type" = "model_type",
                    "PFT" = "PFT",
                    "var_type" = "var_type",
                    "trait" = "trait",
                    "mean" = "Mean",
                    "se_mean" = "se_mean",
                    "sd" = "SD",
                    "2.5%" = "q025",
                    "50%" = "q500",
                    "97.5%" = "q975",
                    "n_eff" = "n_eff",
                    "Rhat" = "Rhat")
    setnames(result, names(names_dict), names_dict)
    result <- result[, names_dict, with=FALSE]

    # Process the PFT name
    result <- result %>%
        separate(PFT, into=c("Biome", "Function"), sep="_", 
                 extra="merge", remove=FALSE) %>%
        separate(Function, into = c("growth_form", "ps_type", 
                                    "leaf_type", "phenology"),
                 sep = "_", extra = "drop", remove=FALSE) %>%
        setDT()
    result[is.na(ps_type), 
        c("ps_type", "leaf_type", "phenology") := Function]


    return(result)
}
