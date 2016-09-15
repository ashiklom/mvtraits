#' @export
buildModel <- function(model_type, dat, model_name = "testmodel",
                       custom_inputs = list()) {
    if(model_type == "uni") {
        stop("uni not supported yet")
        #modelbuild <- buildModel_uni(dat, custom_inputs)
    } else if (model_type == "multi") {
        modelbuild <- buildModel_multi(dat, custom_inputs)
    } else if (model_type == "hier") {
        stop("hier not supported yet")
        #modelbuild <- buildModel_hier(dat, custom_inputs)
    }
    return(modelbuild)
}

