summary_stats <- function(x) {
    out <- c(mean(x), sd(x), quantile(x, c(.025, 0.5, 0.975)))
    names(out) <- c("Mean", "SD", "q025", "q500", "q975")
    return(out)
}

summarizeSigma <- function(samples, pft_name){
    apply(samples, 2, summary_stats) %>%
    t %>%
    "rownames<-"(gsub("Sigma_.*\\.(log\\..*\\.log\\..*)", 
                      "\\1", rownames(.))) %>%
    .[apply(traits.combine, 2, paste, collapse = "."),] %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames("rn", "trait") %>%
    .[, PFT := pft_name]
}

summarizeSigmaPFT <- function(dat, Cor = FALSE) {
    dat.list <- list()
    for (pft_name in pft.names) {
        dat.pft <- getPFTFromHier(dat, pft_name)
        if (Cor) dat.pft <- covToCor.global(dat.pft)
        dat.summary <- summarizeSigma(dat.pft, pft_name)
        dat.list[[pft_name]] <- dat.summary
    }
    out <- rbindlist(dat.list) %>%
        separate(PFT, into=c("Biome", "Function"), sep="_", 
                 extra="merge", remove=FALSE) %>%
        separate(Function, into = c("growth_form", "ps_type", 
                                    "leaf_type", "phenology"),
                sep = "_", extra = "drop", remove=FALSE) %>%
        setDT()
    out[is.na(ps_type), 
        c("ps_type", "leaf_type", "phenology") := Function]
return(out)
}
