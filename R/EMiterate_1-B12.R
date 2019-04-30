#' single iteration of EM algorithm
#'
#' @param Sig list with each element of list being an NxN sigma matrix
#' @param lMS list with lMS[[1]] = M an [N,N,T,T] array and lMS[[2]] = S[T,N]
#' matrix
#'
#' @return updated values of all inputs
#' @export
#'

EMiterate_1_B12 = function(Sig, lMS){
  J = length(Sig)
  for(j in 1:J)  Sig[[j]] = EMcritical(j, Sig, lMS)
  param = sig2param(Sig)
  lMS = sigexRun(param)
  return(list(Sig, lMS))
}
