#' @export
summarizeStan <- function(filename) {
    library(rstan)
    library(data.table)

    # Parse model string
    model_string <- ".*[/](uni|multi|hier)_?([[:digit:]]{0,2})\\.Rdata"
    model_type <- gsub(model_string, "\\1", filename)
    model_pft <- as.numeric(gsub(model_string, "\\2", filename))
    if (is.na(model_pft)) model_pft <- 0

    load(filename)
    outsum <- summary(out)$summary
    result <- as.data.table(outsum, keep.rownames = TRUE)
    rm(out); gc()
    result[, c("model_type", "model_pft") := list(model_type, model_pft)]

    # Parse parameter type
    rnstring <- "^(.*)\\[([[:digit:]]),?([[:digit:]])?,?([[:digit:]])?\\]"
    result[, var_type := gsub(rnstring, "\\1", rn)]

    getnum <- function(num, string) {
        num <- num + 1
        as.numeric(gsub(rnstring, paste0("\\", num), string))
    }

    cols <- c("trait", "var_pft")
    sg <- c("Sigma", "Sigma_global")
    og <- c("Omega", "Omega_global")
    sog <- c(sg, og)
    sop <- c("Sigma_pft", "Omega_pft")

    result[,n1 := getnum(1,rn)][,n2 := getnum(2,rn)][,n3 := getnum(3,rn)]
    result[var_type %in% c("mu", "mu_global"), 
           c(cols) := list(traits_nolog[n1], 0)]

    # Remove duplicates of symmetric matrices
    result <- result[!(var_type %in% sg & n1 > n2)]
    result <- result[!(var_type %in% og & n1 >= n2)]

    ptrait <- function(n1, n2) paste(traits_nolog[n1], traits_nolog[n2], sep=".")

    # Assign Sigma and Omega trait names
    result[var_type %in% sog, c(cols) := list(ptrait(n1, n2), 0)]

    # Hierarchical indices
    result[var_type == "mu_pft", c(cols) := list(traits_nolog[n2], n1)]
    result <- result[!(var_type=="Sigma_pft" & n2 > n3)]
    result <- result[!(var_type=="Omega_pft" & n2 >= n3)]
    result[var_type %in% sop, c(cols) := list(ptrait(n2, n3), n1)]
    
    result <- result[!is.na(trait)]

    names_dict <- c("model_type" = "model_type",
                    "model_pft" = "model_pft",
                    "var_type" = "var_type",
                    "var_pft" = "var_pft",
                    "trait" = "trait",
                    "mean" = "Mean",
                    "se_mean" = "se_mean",
                    "sd" = "SD",
                    "2.5%" = "q025",
                    "50%" = "q500",
                    "97.5%" = "q975")

    setnames(result, names(names_dict), names_dict)
    result <- result[, names_dict, with=FALSE]
    return(result)
}
