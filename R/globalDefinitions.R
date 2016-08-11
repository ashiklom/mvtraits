########################################
# Useful variable definitions
########################################

traits <- c("log.LL", "log.LMA", "log.Nmass", "log.Pmass", "log.Rdmass")
trait.combine <- combn(traits, 2)

try.pfts <- readRDS(system.file("extdata", "try.pft.table.rds", package = "mvtraits"))
pft.factor <- try.pfts$pft.factor
pft.names <- levels(pft.factor)
pft.names[pft.names == "CAM"] <- "arid_CAM"
npft <- length(pft.names)
