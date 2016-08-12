########################################
# Helper functions
########################################

nan2na <- function(x){
  x[is.nan(x)] <- NA
  return(x)
}

is.error <- function(x) class(x) == "try-error"

replace.with.global <- function(x, global){
    nas <- which(is.na(x) | is.nan(x))
    x[nas] <- global[nas]
    return(x)
}

convert.rownames <- function(dat){
    pattern <- "(.*)\\.na\\..*"
    dat$Model <- gsub(pattern, "\\1", rownames(dat))
    rownames(dat) <- NULL
    return(dat)
}

mypng <- function(path, width=800, height=width, ...){
    png(path, width=width, height = height, ...)
}

