customJAGS <- function(model, data, inits, n.chains, 
                       variable.names){
    message("Running model...")
    result <- autorun.jags(model = model,
                           monitor = variable.names,
                           data = data,
                           n.chains = n.chains,
                           inits = inits,
                           max.time = "48h",
                           method = "rjparallel",
                           raftery.options = FALSE)
    return(result)
}

