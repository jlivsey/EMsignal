#' convert param list to Sigma matrix list
#'
#' @param param parameter list in sigex form
#'
#' @return list of NxN Sigma matricies
#' @export
#'

param2sig = function(param){
  J = length(param[[1]])
  Sig = param[[1]] # initialize list
  for(j in 1:J){
      Sig[[j]] = param[[1]][[j]] %*%
                 diag(exp(param[[2]][[j]])) %*%
                 t(param[[1]][[j]])
  }
  return(Sig)
}
