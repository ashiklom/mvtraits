loadTRYData <- function(try_path="data/try.data.rds"){
    try_raw <- readRDS(try_path)
    try_full <- try_raw[,lapply(.SD, nan2na)]
    try_data <-  try_full[,c(traits,"pft.factor"), with = F]
    setnames(try_data, c(traits, "pft"))
    try_data <- try_data[!is.na(try_full$pft)]
    return(try_data)
}
