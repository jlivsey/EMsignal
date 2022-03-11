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

  # put together differencing operators
  diff.full = sigex.delta(mdl = mdl, omits = NA) # full difference operator
  diff.full = round(diff.full) # this is not needed in updated sigex package
  d.full = length(diff.full)
  # Build differencing matrix
  delta0pad = c(diff.full, rep(0, TT-d.full+1))
  D = suppressWarnings(matrix(delta0pad, nrow = TT-d.full,
                              ncol = TT, byrow = TRUE))

  # put together overdifferencing operator coef vectors
  diff.over = list() # over differenced component operators
  for(j in 1:J) diff.over[[j]] = sigex.delta(mdl = mdl, omits = j)

  # form output into paper notation
  M = list()
  M.diff = list()
  for(j in 1:J){

    print(j)

    M[[j]] = block2array(signal[[j]][[2]], N = N, TT = TT)
    # need to define d the full differencing order
    M.diff[[j]] = array(NA, c(N, TT-d.full, N, TT-d.full))
    for(k in (d.full+1):(TT)){
    for(el in (d.full+1):(TT)){
      for(s in 1:d.full){
      for(t in 1:d.full){
        #print(c(k, el, s, t))

        M.diff[[j]][,k-d.full,, el-d.full] = diff.full[s] * diff.full[t] *
                                              M[[j]][,k-s,,el-t]
      }}
    }}
  }


  S = list()
  for(j in 1:J){
    S[[j]] = extract[[j]][[1]]
    S[[j]] = D %*% S[[j]]
  }

  return(list(M, S))
}
