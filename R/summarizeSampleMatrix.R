## Correlation and Covariance Matrices

summarizeSampleMatrix <- function(cov.all.samples, dims, dim.names){
  # Depends on global variables: 
  #     (00.common.R) trait, trait.combine

  # Calculate summary statistics across samples
  cov.all.list <- list(Mean = apply(cov.all.samples, dims, mean),
                       SD = apply(cov.all.samples, dims, sd),
                       q025 = apply(cov.all.samples, dims, quantile, 0.025),
                       q500 = apply(cov.all.samples, dims, quantile, 0.500),
                       q975 = apply(cov.all.samples, dims, quantile, 0.975)
  ) %>% lapply("dimnames<-",  dim.names)
  
  getcov <- function(trait, cov.all) cov.all[, trait[1], trait[2]]
  columnize <- function(mat.wide){
    cov.mat <- apply(trait.combine, 2, getcov, mat.wide) %>% as.data.frame
    cov.mat.names <- apply(trait.combine, 2, paste, collapse="_") %>%
      gsub("log.", "", .)
    colnames(cov.mat) <- cov.mat.names
    return(cov.mat)
  }
  
  cov.dat <- lapply(cov.all.list, columnize) %>% 
    do.call(cbind, .) %>%
    add_rownames(var = "PFT") %>%
    gather(Stat.Trait, Value, -PFT) %>%
    separate(Stat.Trait, into=c("Stat", "Trait"), sep="\\.") %>%
    spread(Stat, Value) %>%
    separate(PFT, into=c("Biome", "Function"), sep="_", 
             extra="merge", remove=FALSE)
  
  return(cov.dat)
}
