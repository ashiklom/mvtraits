#' @export
buildModel <- function(model_type, dat, model_name = "testmodel",
                       custom_inputs = list()) {
    if(model_type == "uni") {
        modelbuild <- buildModel_uni(dat, custom_inputs)
    } else if (model_type == "multi") {
        modelbuild <- buildModel_multi(dat, custom_inputs)
    } else if (model_type == "hier") {
        modelbuild <- buildModel_hier(dat, custom_inputs)
    }
    return(modelbuild)
}

