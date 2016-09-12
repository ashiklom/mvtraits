# Functions for easy operations on multidimensional arrays

#' @export
array2DApply <- function(arr, func, ...){
    arr.dim <- dim(arr)
    out <- array(NA, arr.dim)
    for(i in 1:arr.dim[1]) out[i,,] <- func(arr[i,,], ...)
    return(out)
}

#' @export
array3Dapply <- function(arr, func, ...){
    out <- array(NA, dim(arr))
    for(i in 1:dim(out)[1]){
        for(j in 1:dim(out)[2]){
            out[i,j,,] <- func(arr[i,j,,], ...)
        }
    }
    return(out)
}

