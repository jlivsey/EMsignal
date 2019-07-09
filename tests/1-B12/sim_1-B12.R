library(sigex)
library(mvtnorm)
library(EMsigex)

# ---- Simulate Data ----------------------------------------------------------
N = 3
T <- TT <- 300
t = 1:T
Phi=diag(N)
Sig=diag(N); Sig[1,2] <- Sig[2,1] <- .75

s1 = gen_trendComp(T, Phi, Sig)
s2 = gen_seasComp(T, Phi, diag(N))
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
data = demean(data)

plot(ts(data), main="simulated series")

# ---- plotting and checks ----------------------------------------------------
# plot(ts(data))
# acf(data, lag.max = 100)
# acf(diff(s1))
# acf(diff(s2, 12))
# acf(s0)
par(mfrow=c(dim(data)[2],1))
spec.ar(ts(data[,1]))
spec.ar(ts(data[,2]))
spec.ar(ts(data[,3]))

# ---- Modeling ---------------------------------------------------------------
transform = "none"
x <- t(data)
N <- dim(x)[1]
T <- dim(x)[2]

# ---- Model ------------------------------------------------------------------
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"seasonal", rep(1,12))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl,data,0)



# Set default parameters
par.default <- sigex.default(mdl,data)[[1]]
flag.default <- sigex.default(mdl,data)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)

# Set param to TRUE values
param = par.default

# # ---- MOM estimates for param ------------------------------------------------
mdl.mom <- mdl
par.mom <- sigex.momfit(data,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,flag.default,mdl.mom)
#resid.mom <- sigex.resid(psi.mom,mdl.mom,data)

thresh <- -1.66

if(N > 1) {
  reduced.mom <- sigex.reduce(data,par.mom,flag.default,mdl.mom,thresh,FALSE)
  mdl.mom <- reduced.mom[[1]]
  par.mom <- reduced.mom[[2]]
  flag.mom <- sigex.default(mdl.mom,data)[[2]]
  psi.mom <- sigex.par2psi(par.mom,flag.mom,mdl.mom)
  #resid.mom <- sigex.resid(psi.mom,mdl.mom,data)
}

# bundle for default span
analysis.mom <- sigex.bundle(data,transform,mdl.mom,psi.mom)

## Rough: reduced MOM model
#data <- analysis.mom[[1]]
#mdl <- analysis.mom[[3]]
psi <- analysis.mom[[4]]

param.mom <- sigex.psi2par(psi,mdl,data)

# ---- DMA for model ----------------------------------------------------------
signal.trendann <- sigex.signal(data,param,mdl,1)
signal.seas     <- sigex.signal(data,param,mdl,2)
signal.irr      <- sigex.signal(data,param,mdl,3)

extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
extract.seas     <- sigex.extract(data,signal.seas,mdl,param)
extract.irr      <- sigex.extract(data,signal.irr,mdl,param)


# ---- Plots ------------------------------------------------------------------
subseries <- 2
xss = data[, subseries]
s1.hat  = extract.trendann[[1]][, subseries]
s2.hat  = extract.seas[[1]][, subseries]
s0.hat = extract.irr[[1]][, subseries]
{
  op = par(mfrow=c(3,1), mar=c(2,3,2,1))
  plot(as.numeric(xss), type="l")
  lines(s1.hat, col="tomato")
  plot(s2.hat, type="l", col="seagreen"); abline(h=0, lty="dotted")
  abline(v=seq(1,TT,12), lty="dashed")
  plot(s0.hat, type="l", col="navyblue"); abline(h=0, lty="dotted")
  par(op)
}

# --- Filter weights ----------------------------------------------------------
FF = signal.trendann[[1]]
FF = block2array(FF, N, T)
fw = FF[1, T/2, 1, ]
plot(fw, type="l")
abline(v=seq(T/2, T, 12), lty="dotted")
sum(fw)

