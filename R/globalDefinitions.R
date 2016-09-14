########################################
# Useful variable definitions
########################################

#' @export
traits <- c("log.LL", "log.LMA", "log.Nmass", "log.Pmass", "log.Rdmass")

#' @export
traits_nolog <- gsub("log\\.", "", traits)

#' @export
trait.combine <- combn(traits, 2)

#' @export
try.pfts <- readRDS(system.file("extdata", "try.pft.table.rds", package = "mvtraits"))

#' @export
pft.factor <- try.pfts$pft.factor

#' @export
pft.names <- levels(pft.factor)

#' @export
pft.names[pft.names == "CAM"] <- "arid_CAM"

#' @export
npft <- length(pft.names)
