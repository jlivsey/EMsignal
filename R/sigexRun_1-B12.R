#' Runs sigex for EM algorithm
#'
#' @param param parameter list in sigex form
#' @param data  mts object (sample size) x (dim)
#' @param mdl   mdl object in sigex form
#'
#' @return list with M (error matricies) and S (signal estimates)
#' @export
#'

sigexRun_1_B12 = function(param, data, mdl){

  TT <- dim(data)[1]
  N  <- dim(data)[2]

  # Run sigex with updated param
  signal.trendann <- sigex.signal(data,param,mdl,1)
  signal.seas     <- sigex.signal(data,param,mdl,2)
  signal.irr      <- sigex.signal(data,param,mdl,3)

  extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
  extract.seas     <- sigex.extract(data,signal.seas,mdl,param)
  extract.irr      <- sigex.extract(data,signal.irr,mdl,param)

  # form output into paper notation
  M1 = block2array(signal.trendann[[2]], N = N, TT = TT)
  M2 = block2array(signal.seas[[2]],     N = N, TT = TT)
  M3 = block2array(signal.irr[[2]],      N = N, TT = TT)
  M = list(M1, M2, M3)

  S1 = extract.trendann[[1]]
  S1d = diff(S1, 12)
  S2 = extract.seas[[1]]
  S2d = diff(S2, 12)
  S3 = extract.irr[[1]]
  S3d = diff(S3, 12)
  S = list(S1d, S2d, S3d)

  return(list(M, S))
}
