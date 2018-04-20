#' Convert list of Sigma matricies to sigex param list
#'
#' @param Sig list with each element of list being an NxN sigma matrix
#'
#' @return updated param list
#' @export
#'

sig2param = function(Sig){
  J = length(Sig)
  for(j in 1:J){
    GCD = getGCD(Sig[[j]], 2)
    param[[1]][[j]] = GCD[[1]]
    param[[2]][[j]] = log(GCD[[2]])
  }
  return(param)
}
