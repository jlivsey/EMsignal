
#' Generate multivariate trend component
#'
#' generates a component that is first difference stationary
#'
#' @param n length of series
#' @param Phi autoregressive parameter
#' @param Sig covariance martrix white noise component
#' @param burn burn in (defaults to 1000)
#'
#' @return Ndim x n matrix of observations
#' @export
#'

gen_trendComp = function(n, Phi, Sig, burn=1000){
  N = n+burn
  Ndim = dim(Sig)[1]

  if(Ndim==1){ # handle univariate case first
    w = rnorm(n = N, mean = 0, sd = as.numeric(sqrt(Sig)))
    s = w[1]
    for(i in 2:(N)){
      new.s = Phi * s[i-1] - w[i]
      s = c(s, new.s)
    }
    return(s[(burn+1):N])
  } # end of univariate if() statement

  w = rmvnorm(n = N+1, mean = rep(0,Ndim), sigma = Sig)
  s = w[1:2, ] # Need at least 2 here or don't get matrix class object
  for(i in 3:(N+1)){
    new.s = Phi %*% s[i-1, ] - w[i, ]
    s = rbind(s, t(new.s))
  }
  return(s[(burn+1):N, ])
}
