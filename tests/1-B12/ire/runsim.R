# library(sigex)
library(mvtnorm)
# library(EMsigex)
library(parallel)

# ---- Load data ----
load("dataList.Rdata")

# ---- Simulation function -----------------------------------------------------

sim <- function(data.ts){

  tryCatch(
    expr = {

# ---- Modeling ---------------------------------------------------------------
transform = "none"
T <- dim(data.ts)[1]
TT <- T
N <- dim(data.ts)[2]

# ---- Model ------------------------------------------------------------------
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"seasonal", rep(1,12))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl, data.ts, 0)

# Set default parameters
constraint <- NULL
par.default <- sigex.default(mdl, data.ts, constraint)
psi.default <- sigex.par2psi(par.default, mdl)


# Set param to TRUE values
param = par.default

# # ---- MOM estimates for param ------------------------------------------------
mdl.mom <- mdl
par.mom <- sigex.momfit(data.ts, par.default, mdl.mom)
psi.mom <- sigex.par2psi(par.mom, mdl.mom)

thresh <- -1.66
if(N > 1) {
  reduced.mom <- sigex.reduce(data.ts = data.ts,
                              param = par.mom,
                              mdl = mdl.mom,
                              thresh = thresh,
                              modelflag = FALSE)
  mdl.mom <- reduced.mom[[1]]
  par.mom <- reduced.mom[[2]]
  psi.mom <- sigex.par2psi(par.mom,mdl.mom)
}

# bundle and extract psi/param
analysis.mom <- sigex.bundle(data.ts ,transform, mdl.mom, psi.mom)
psi          <- analysis.mom[[4]]
par.mom      <- sigex.psi2par(psi, mdl, data.ts)

# ---- MLE --------------------------------------------------------------------

## load up the MOM model
data.ts <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi.mom <- analysis.mom[[4]]
par.mom <- sigex.psi2par(psi.mom,mdl,data.ts)

#  Initialize with MOM estimates
constraint <- NULL
psi.mle <- sigex.par2psi(par.mom, mdl)

## run fitting: can be commented out, this takes a while
fit.mle <- sigex.mlefit(data.ts,
                        par.mom,
                        constraint,
                        mdl,
                        method = "bfgs",
                        debug  = FALSE)

# set psi from mle estimation
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian

# bundle for default span
analysis.mle <- sigex.bundle(data.ts, transform, mdl, psi.mle)
par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)

# ---- signal extraction from fit - Direct Matrix Approach ---------------------

# which params to use default/mom/mle
param <- par.mom

signal.trendann <- sigex.signal(data.ts, param, mdl, 1)
signal.seas     <- sigex.signal(data.ts, param, mdl, 2)
signal.irr      <- sigex.signal(data.ts, param, mdl, 3)

extract.trendann <- sigex.extract(data.ts, signal.trendann, mdl, param)
extract.seas     <- sigex.extract(data.ts, signal.seas,     mdl, param)
extract.irr      <- sigex.extract(data.ts, signal.irr,      mdl, param)



# ------------------------------------------------------------------------------
# ---- Start EM evaluation -----------------------------------------------------
# ------------------------------------------------------------------------------


d = 12
Gam1 = toeplitz(ARMAauto(ma = rep(1,11), ar = NULL, lag.max = (TT-1-d)))
Gam2 = toeplitz(ARMAauto(ma = -1, ar = NULL, lag.max = (TT-1-d)))
Gam3 = toeplitz(ARMAauto(ma = c(rep(0,11), -1), ar = NULL, lag.max = (TT-1-d)))
Gam  = list(Gam1, Gam2, Gam3)
invGam = lapply(Gam, solve)

# ---- Likelihood at the TRUE values ------------------------------------------

Sig1 <- diag(N)
Sig1[1,2] <- Sig1[2,1] <- .75
Sig.true <- list(Sig1, diag(N), diag(N))
lik.true <- sig2lik(Sig.true, mdl, data.ts)

# ---- Likelihood at the MOM  ------------------------------------------

param = par.mom
Sig.mom = param2sig(param)
lik.mom <- sig2lik(Sig.mom, mdl, data.ts)

# ---- Likelihood at the MLE  ------------------------------------------

param = par.mle
Sig.mle = param2sig(param)
lik.mle <- sig2lik(Sig.mle, mdl, data.ts)

# ---- Initialize values for first iteration ----------------------------------

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

# ---- Run EM ----------------------------------------


iters <- 50
Nc <- length(unlist(Sig.mom))
Sig.save <- matrix(NA, nrow = iters+2, ncol= Nc+1)
Sig.save[1, ] <- c(unlist(Sig.true), lik.true)
Sig.save[2, ] <- c(unlist(Sig.mom), lik.mom)
for(i in 1:iters) {
    if(i==1) Sig <- Sig.mom
  out = EMiterate_1_B12(Sig, lMS, data.ts, mdl)
  Sig = out[[1]]
  lMS = out[[2]]
  lik <- sig2lik(Sig, mdl, data.ts)
  Sig.save[i+2, ] <- c(unlist(Sig), lik)
  # print(Sys.time())
  # print(i)
  # print(lik)
  # print("------------------------------------")
}


return(list(lik.true = lik.true,
            lik.mom  = lik.mom,
            lik.mle  = lik.mle,
            psi.mle  = psi.mle,
            psi.mom  = psi.mom,
            Sig.save = Sig.save))

  }, # END of expr = {
  error = function(e){
    message('Caught an error!')
    print(e)
  },
  warning = function(w){
    message('Caught an warning!')
    print(w)
  }
)


}

c <- detectCores()
out <- mclapply(dataList, sim, mc.cores = c)
# save(out, file = 'out.Rdata')







