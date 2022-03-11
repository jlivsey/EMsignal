#' Calculate updated Sigma values at EM critical point
#'
#' @param j component of mdl
#' @param lMS list with lMS[[1]] = M an [N,N,T,T] array and lMS[[2]] = S[T,N]
#' matrix
#'
#' @return NxN updated Sigma matrix
#' @export
#'

EMcritical = function(j, Sig, lMS, mdl){
  M = lMS[[1]]
  S = lMS[[2]]
  #d = length(mdl[[3]][[j]])
  d <- 12

  # put together overdifferencing operator coef vectors
  diff.over = sigex.delta(mdl = mdl, omits = j)

  for(k in (d+1):TT){
    for(el in (d+1):TT){
      Mj = M[[j]][,k,,el] - M[[j]][,k-d,,el] -
           M[[j]][,k,,el-d] + M[[j]][,k-d,,el-d]
    }
  }

  outSig = matrix(0, N, N)
  for(k in (d+1):TT){
    for(el in (d+1):TT){
      D = invGam[[j]][k-d, el-d] # * invSig[[1]]
      new.term =  D * (Mj + (S[[j]][k-d, ] %*% t(S[[j]][el-d, ])))
      outSig = outSig + new.term
    }}
  outSig = outSig / (TT-d)
  return(outSig)
  #return(( outSig + t(outSig) )/2)
}
