source("sim1-B12.R")

Gam1 = toeplitz(ARMAacvf(ma=rep(-1,11), lag.max = (TT-1)))
Gam2 = toeplitz(ARMAacvf(ma=1, lag.max = (TT-1)))
Gam3 = toeplitz(ARMAacvf(ma=c(rep(0,11), 1), lag.max = (TT-1)))
Gam  = list(Gam1, Gam2, Gam3)
invGam = lapply(Gam, solve)

param = par.default
param <- sigex.psi2par(psi,mdl,data)

# ---- Initialize values for first iteration ----------------------------------

Sig = param2sig(param)

M1 = block2array(signal.trendann[[2]], N = N, TT = TT)
M2 = block2array(signal.seas[[2]],     N = N, TT = TT)
M3 = block2array(signal.irr[[2]],      N = N, TT = TT)
M = list(M1, M2, M3)

S1 = extract.trendann[[1]]
S2 = extract.seas[[1]]
S3 = extract.irr[[1]]
S = list(S1, S2, S3)

lMS = list(M, S)

# -----------------------------------------------------------------------------

out = EMiterate(Sig, lMS)
(Sig = out[[1]])
lMS = out[[2]]

out = EMiterate(Sig, lMS)
(Sig = out[[1]])
lMS = out[[2]]

