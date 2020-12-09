#' Generate Data from Trend + Seasonal + Irregular model
#' Assumes (1 - B) reduces trend to stationary and (1-B)^{12} reduces
#' the seasonal component to stationary
#'
#' @param N dimensionality
#' @param TT sample size
#'
#' @return multivariate time series of sigex form
#' @export
#'

gen_data <- function(N, TT){

  start.date <- c(1990, 1) # generic
  period     <- 12
  transform  <- "none"
  aggregate  <- FALSE

  Phi <- diag(N)
  Sig <- diag(N); Sig[1,2] <- Sig[2,1] <- .75

  # ---- generate structural components
  s1 = gen_trendComp(TT, Phi, Sig)
  s2 = gen_seasComp(TT, Phi, diag(N))
  s0 = mvtnorm::rmvnorm(n = TT, mean = rep(0,N), sigma = diag(N))

  data = s1+s2+s0
  data = demean(data)

  # ---- create dim names for sigex.load ----
  dimNames <- rep(NA, N)
  for(i in 1:N) dimNames[i] <- paste0("dim", i)

  # ---- Load ----
  dataALL.ts <-
    sigex::sigex.load(data = data,
                      start.date = start.date,
                      period = period,
                      epithets = dimNames,
                      plot = FALSE)

  # ---- Prep for sigex ----
  begin.date <- start.date
  end.date <- end(dataALL.ts)
  range <- list(begin.date,end.date)
  data.ts <-
    sigex::sigex.prep(data.ts = dataALL.ts,
                      transform = transform,
                      aggregate = aggregate,
                      subseries = 1:N,
                      range = range,
                      plot = FALSE)

  return(data.ts)

}
