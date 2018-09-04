

#' Demean vector time series
#'
#' takes a sample of size (T) x dim (N) matrix (or its transpose) and removes the
#' sample means of each column (row).
#'
#' @param x matrix vector valued time series (assumes dim < sample size)
#'
#' @return same dim matrix as supplied with sample mean removed
#' @export
#'

demean = function(x){
  if(dim(x)[1] > dim(x)[2]){
    x = t(x)
    return(t(x - rowMeans(x)))
  } else{
    return(x - rowMeans(x))
  }
}



