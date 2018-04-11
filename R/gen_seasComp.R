
#' Generate multivariate seasonal component
#'
#' generates a component that is seasonally (12) difference stationary
#'
#' @param n length of series
#' @param Phi seasonal autoregressive parameter
#' @param Sig covariance martrix white noise component
#' @param burn burn in (defaults to 1000)
#' @param seasonl.period periodicity of component (defaults to 12)
#'
#' @return Ndim x n matrix of observations
#' @export
#'

gen_seasComp = function(n, Phi, Sig, burn=10^3, seasonal.period=12){
  N = n+burn
  Ndim = dim(Sig)[1]
  w = rmvnorm(n = N, mean = rep(0,Ndim), sigma = Sig)
  s = w[1:seasonal.period, ]
  print(s)
  # for(i in (seasonal.period+1):N){
  #     new.s = Phi %*% s[i-seasonal.period, ] - w[i, ]
  #     s = rbind(s, t(new.s))
  # }
  for(i in (seasonal.period+1):N){
    new.s = -1*colSums(s[(i-seasonal.period):(i-1), ]) + w[i,]
    s = rbind(s, t(new.s))
  }
  return(s[(burn+1):N, ])
}
