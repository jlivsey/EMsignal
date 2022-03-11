#' Perform multivariate signal extraction with *sigex* via the EM algoritm
#'
#' @param data multivariate time series as a matrix (sampleSize x dim)
#' @param mdl sigex model object
#' @param transform transform for data. Default "none"
#'
#' @return list of parameter estimates
#' @export
#'

EMsigex <- function(data, mdl, transform="none"){

  J = length(mdl[[2]]) # total number of components in model
  d = unlist(lapply(mdl[[3]], length)) # diff order of each component

  # put together overdifferencing operator coef vectors
  diff.over = list()
  for(j in 1:J) diff.over[[j]] = sigex.delta(mdl = mdl, omits = j)

  x <- t(data)
  N <- dim(x)[1]
  TT <- dim(x)[2]

  # Set default parameters
  par.default <- sigex.default(mdl,data)[[1]]
  flag.default <- sigex.default(mdl,data)[[2]]
  psi.default <- sigex.par2psi(par.default,flag.default,mdl)

  # Set param to TRUE values
  param = par.default

  # # ---- MOM estimates for param --------------------------------------------
  mdl.mom <- mdl
  par.mom <- sigex.momfit(data,par.default,mdl.mom)
  psi.mom <- sigex.par2psi(par.mom,flag.default,mdl.mom)
  #resid.mom <- sigex.resid(psi.mom,mdl.mom,data)

  thresh <- -1.66

  if(N > 1) {
    reduced.mom <- sigex.reduce(data,par.mom,flag.default,mdl.mom,thresh,FALSE)
    mdl.mom <- reduced.mom[[1]]
    par.mom <- reduced.mom[[2]]
    flag.mom <- sigex.default(mdl.mom,data)[[2]]
    psi.mom <- sigex.par2psi(par.mom,flag.mom,mdl.mom)
    #resid.mom <- sigex.resid(psi.mom,mdl.mom,data)
  }

  # bundle for default span
  analysis.mom <- sigex.bundle(data,transform,mdl.mom,psi.mom)

  ## Rough: reduced MOM model
  #data <- analysis.mom[[1]]
  #mdl <- analysis.mom[[3]]
  psi <- analysis.mom[[4]]

  param.mom <- sigex.psi2par(psi,mdl,data)

  # ---- Initialize values for first iteration ---------------------------------

  param = param.mom

  lMS = sigexRun(param = param, data = data, mdl = mdl)

  # ---- Run EMiterate ----

  for(i in 1:10) {
    out = EMiterate(Sig, lMS); (Sig = out[[1]]); lMS = out[[2]]
    A = labsdiff(Sig, Sig.mle)
    print(A)
    print(sum(unlist(lapply(A, sum))))
    print("------------------------------------")
  }
}
