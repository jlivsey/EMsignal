#' convert param list to correlation matrix list
#'
#' @param param parameter list in sigex form
#'
#' @return list of NxN correlation matricies
#' @export
#'

param2cor = function(param){
  J = length(param[[1]])
  Sig = param[[1]] # initialize list
  for(j in 1:J){
      Sig[[j]] = param[[1]][[j]] %*%
                 diag(exp(param[[2]][[j]])) %*%
                 t(param[[1]][[j]])
  }
  R = lapply(Sig, cov2cor)
  return(R)
}
