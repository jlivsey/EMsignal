# source("sim1-B12.R")

d = 12
Gam1 = toeplitz(ARMA2acf(ma=rep(1,11), lag.max = (TT-1-d)))
Gam2 = toeplitz(ARMA2acf(ma=-1, lag.max = (TT-1-d)))
Gam3 = toeplitz(ARMA2acf(ma=c(rep(0,11), -1), lag.max = (TT-1-d)))
Gam  = list(Gam1, Gam2, Gam3)
invGam = lapply(Gam, solve)

param = param.mom
#param <- sigex.psi2par(psi,mdl,data)

# ---- Initialize values for first iteration ----------------------------------

Sig = param2sig(param)

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

lMS = list(M, S)

# -----------------------------------------------------------------------------

Sig.mle = param2sig(param.mle)
Sig.mle = param2sig(param.mom)
Sig.mle = param2sig(param)


for(i in 1:10) {

  i <- 1

  out = EMiterate(Sig, lMS); (Sig = out[[1]]); lMS = out[[2]]
  A = labsdiff(Sig, Sig.mle)
  print(A)
  print(sum(unlist(lapply(A, sum))))
  print("------------------------------------")
}
