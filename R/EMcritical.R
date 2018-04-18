Gam1 = toeplitz(ARMAacvf(ma=rep(-1,11), lag.max = (TT-1)))
Gam2 = toeplitz(ARMAacvf(ma=1, lag.max = (TT-1)))
Gam3 = toeplitz(ARMAacvf(ma=c(rep(0,11), 1), lag.max = (TT-1)))
Gam = list(Gam1, Gam2, Gam3)

Sig1 <- Sig2 <- Sig3 <- diag(N)
Sig = list(Sig1, Sig2, Sig3)

M1 = block2array(signal.trendann[[2]], N = N, TT = TT)
M2 = block2array(signal.seas[[2]],     N = N, TT = TT)
M3 = block2array(signal.irr[[2]],      N = N, TT = TT)
M = list(M1, M2, M3)

S1 = extract.trendann[[1]]
S2 = extract.seas[[1]]
S3 = extract.irr[[1]]
S = list(S1, S2, S3)

EMcritical = function(j){
  outSig = matrix(0, N, N)
  for(k in 1:TT){
  for(el in 1:TT){
    new.term = Gam[[j]][k, el] * Sig[[1]] %*%
               (M[[j]][,,k,el] + S[[j]][k, ] %*% t(S[[j]][el, ]))
    outSig = outSig + new.term
  }}
  outSig = outSig / TT
  return(outSig)
}

EMcritical(1)
EMcritical(2)
EMcritical(3)

