#' single iteration of EM algorithm
#'
#' @param Sig list with each element of list being an NxN sigma matrix
#' @param lMS list with lMS[[1]] = M an [N,N,T,T] array and lMS[[2]] = S[T,N]
#' matrix
#'
#' @return updated values of all inputs
#' @export
#'

EMiterate = function(Sig, lMS, data, mdl){
  J = length(Sig)
  for(j in 1:J){
    Sig[[j]] = EMcritical(j, Sig, lMS, mdl)
  }
  param = sig2param(Sig)
  lMS = sigexRun(param, data, mdl)
  return(list(Sig, lMS))
}
