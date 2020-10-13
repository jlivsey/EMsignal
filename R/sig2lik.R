#' output likelihood value for given Sigma values
#'
#' @param Sig List of Sigma's to output likelihood of
#'
#' @return likelihood value
#' @export
#'

sig2lik <- function(Sig, mdl = mdl, data = data){
  # Put in param form
  param <- sig2param(Sig, mdl, data)
  # Put in psi form
  psi <- sigex.par2psi(param, mdl)
  # likelihood
  lik = sigex.lik(psi, mdl, data)

  return(as.numeric(lik))
}
