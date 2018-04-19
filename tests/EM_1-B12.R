Gam1 = toeplitz(ARMAacvf(ma=rep(-1,11), lag.max = (TT-1)))
Gam2 = toeplitz(ARMAacvf(ma=1, lag.max = (TT-1)))
Gam3 = toeplitz(ARMAacvf(ma=c(rep(0,11), 1), lag.max = (TT-1)))
invGam = list(solve(Gam1), solve(Gam2), solve(Gam3))

# ---- Initialize values for first iteration ----------------------------------

# True Sig
Sig1 <- Sig2 <- Sig3 <- diag(N)
Sig = list(Sig1, Sig2, Sig3)
invSig = list(solve(Sig1), solve(Sig2), solve(Sig3))

# Piece together Sig from MOM estimates
Sig = list(param[[1]][[1]] %*% diag(param[[2]][[1]]) %*% t(param[[1]][[1]]),
           param[[1]][[2]] %*% diag(param[[2]][[2]]) %*% t(param[[1]][[2]]),
           param[[1]][[3]] %*% diag(param[[2]][[3]]) %*% t(param[[1]][[3]]))
invSig = lapply(Sig, solve)

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
Sig = out[[1]]
lMS = out[[2]]


