#' output likelihood value for given Sigma values
#'
#' @param Sig List of Sigma's to output likelihood of
#'
#' @return likelihood value
#' @export
#'

sig2lik <- function(Sig){
  # Put in param form
  param <- sig2param(Sig)
  # Put in psi form
  psi <- sigex.par2psi(param, flag.mom, mdl)

  # likelihood
  lik = sigex.lik(psi, mdl, data)

  return(as.numeric(lik))
}
