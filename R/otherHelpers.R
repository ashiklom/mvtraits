########################################
# Helper functions
########################################

#' @export
nan2na <- function(x){
  x[is.nan(x)] <- NA
  return(x)
}

#' @export
is.error <- function(x) class(x) == "try-error"

#' @export
replace.with.global <- function(x, global){
  nas <- which(is.na(x) | is.nan(x))
  x[nas] <- global[nas]
  return(x)
}

#' @export
convert.rownames <- function(dat){
  pattern <- "(.*)\\.na\\..*"
  dat$Model <- gsub(pattern, "\\1", rownames(dat))
  rownames(dat) <- NULL
  return(dat)
}

#' @export
graphic <- function(path, width=800, height=width, ...){
  fname <- paste0(path, ".tiff")
  tiff(fname, width=width, height = height, ...)
}

#' @export
getPFTFromHier <- function(hier_dat, pft_name){
  hier_dat[, grep(pft_name, colnames(hier_dat))]
}
