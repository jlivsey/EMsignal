#' Runs sigex for EM algorithm
#'
#' @param param parameter list in sigex form
#' @param data data matrix (dim x sample_size)
#' @param mdl sigex mdl
#'
#' @return list with M (error matricies) and S (signal estimates)
#' @export
#'

sigexRun = function(param, data, mdl){

  N = dim(data)[2]
  TT = dim(data)[1]

  J = length(mdl[[2]]) # total number of components in model
  d = lapply(mdl[[3]], length) # diff order of each component

  signal = rep(list(matrix(-99, N*TT, N*TT)), J) # initalize storage
  for(j in 1:J){
    signal[[j]] = sigex.signal(data, param, mdl, j)
  }

  extract = rep(list(matrix(-99, TT, N)), J) # initalize storage
  for(j in 1:J){
    extract[[j]] = sigex.extract(data, signal[[j]], mdl, param)
  }

  # put together overdifferencing operator coef vectors
  diff.full = sigex.delta(mdl = mdl, omits = NA) # full difference operator
  diff.over = list() # over differenced component operators
  for(j in 1:J) diff.over[[j]] = sigex.delta(mdl = mdl, omits = j)

  # form output into paper notation
  M = list()
  for(j in 1:J){
    M[[j]] = block2array(signal[[j]][[2]], N = N, TT = TT)
  # differenced M double sum goes here !!!
  }

  S = list()
  for(j in 1:J){
    S[[j]] = extract[[j]][[1]]
    S[[j]] = filter(x = S[[j]], filter = diff.over[[j]], sides = 1)
    S[[j]] = na.omit(S[[j]])
    S[[j]] = as.matrix(S[[j]])
  }

  return(list(M, S))
}
