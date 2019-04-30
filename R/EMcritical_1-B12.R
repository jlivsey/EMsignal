#' Calculate updated Sigma values at EM critical point
#'
#' @param j component of mdl
#' @param lMS list with lMS[[1]] = M an [N,N,T,T] array and lMS[[2]] = S[T,N]
#' matrix
#'
#' @return NxN updated Sigma matrix
#' @export
#'

EMcritical_1_B12 = function(j, Sig, lMS){
  M = lMS[[1]]
  S = lMS[[2]]
  d = 12
  invSig = lapply(Sig, solve)
  outSig = matrix(0, N, N)
  for(k in (d+1):TT){
    for(el in (d+1):TT){
      D = invGam[[j]][k-d, el-d] # * invSig[[1]]
      new.term =  D * (M[[j]][,,k,el] - M[[j]][,,k-12,el] - M[[j]][,,k,el-12] +
                         M[[j]][,,k-12,el-12] + (S[[j]][k-d, ] %*% t(S[[j]][el-d, ])))
      outSig = outSig + new.term
    }}
  outSig = outSig / (TT-d)
  return(outSig)
  #return(( outSig + t(outSig) )/2)
}
