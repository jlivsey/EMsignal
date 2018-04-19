#' Calculate updated Sigma values at EM critical point
#'
#' @param j component of mdl
#' @param lMS list with lMS[[1]] = M an [N,N,T,T] array and lMS[[2]] = S[T,N]
#' matrix
#'
#' @return NxN updated Sigma matrix
#' @export
#'

EMcritical = function(j, lMS){
  M = lMS[[1]]
  S = lMS[[2]]
  outSig = matrix(0, N, N)
  for(k in 1:TT){
    for(el in 1:TT){
      D = invGam[[j]][k, el] * invSig[[1]]
      new.term =  D %*% (M[[j]][,,k,el] + S[[j]][k, ] %*% t(S[[j]][el, ]))
      outSig = outSig + new.term
    }}
  outSig = outSig / TT
  return(outSig)
}
