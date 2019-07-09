#' Convert list of Sigma matricies to sigex param list
#'
#' @param Sig list with each element of list being an NxN sigma matrix
#'
#' @return updated param list
#' @export
#'

sig2param = function(Sig){
  J <- length(Sig)
  N <- dim(Sig[[1]])[1]
  for(j in 1:J){
    GCD = getGCD(Sig[[j]], N)
    param[[1]][[j]] = GCD[[1]]
    param[[2]][[j]] = log(GCD[[2]])
  }
  return(param)
}
