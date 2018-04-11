library(sigex)
library(mvtnorm)

# ---- Simulate Data ----------------------------------------------------------
N = 2
T = 200
t = 1:T
Phi=diag(N)
Sig=diag(N)

s1 = gen_trendComp(T, Phi, Sig)
s2 = gen_seasComp(T, Phi, Sig)
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
#plot(ts(data))

#spec.ar(data[,1])
#abline(v=seq(1/12, 1/2, 1/12))

# ---- Modeling ---------------------------------------------------------------
transform = "none"
agg <- FALSE	# set TRUE to aggregate
series <- 1:3
range <- seq(1,dim(data)[1])
x <- t(data)
N <- dim(x)[1]
T <- dim(x)[2]

# ---- Default Model ----------------------------------------------------------

# stochastic effects
delta.trend <- c(1,-1)
def <- c(0,1,0,1)

mdl <- NULL
# trend-cycle (annual cycle)
mdl <- sigex.add(mdl,seq(1,N),"wn",delta.trend,def)
# first atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(2*pi/12),1),def)

# second atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(4*pi/12),1),def)
# third atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(6*pi/12),1),def)
# fourth atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(8*pi/12),1),def)
# firth atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(10*pi/12),1),def)
# pi freq
mdl <- sigex.add(mdl,seq(1,N),"wn", c(1,1), def)
# irregular
mdl <- sigex.add(mdl,seq(1,N),"wn",1,def)
# fixed effects [mean effect]
mdl <- sigex.meaninit(mdl,data,0)

par.default <- sigex.default(mdl,data)[[1]]
flag.default <- sigex.default(mdl,data)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)

# ---- MOM estimation and reduced specification -------------------------------
mdl.mom <- mdl
par.mom <- sigex.momfit(data,par.default,mdl.mom)
psi.mom <- sigex.par2psi(par.mom,flag.default,mdl.mom)
resid.mom <- sigex.resid(psi.mom,mdl.mom,data)

#thresh <- -6.22
#thresh <- -3.92
thresh <- -1.66

if(N > 1) {
  reduced.mom <- sigex.reduce(data,par.mom,flag.default,mdl.mom,thresh,FALSE)
  mdl.mom <- reduced.mom[[1]]
  par.mom <- reduced.mom[[2]]
  flag.mom <- sigex.default(mdl.mom,data)[[2]]
  psi.mom <- sigex.par2psi(par.mom,flag.mom,mdl.mom)
  #resid.mom <- sigex.resid(psi.mom,mdl.mom,data)
}

log(sigex.conditions(data,psi.mom,mdl.mom))


# model checking
sigex.portmanteau(t(resid.mom),48,length(psi.mom))
sigex.gausscheck(t(resid.mom))

# bundle for default span
analysis.mom <- sigex.bundle(data,transform,mdl.mom,psi.mom)

## Rough: reduced MOM model
data <- analysis.mom[[1]]
mdl <- analysis.mom[[3]]
psi <- analysis.mom[[4]]
param <- sigex.psi2par(psi,mdl,data)

# ---- METHOD 1: DIRECT MATRIX APPROACH ---------------------------------------
signal.trendann <- sigex.signal(data,param,mdl,1)
signal.seas.1   <- sigex.signal(data,param,mdl,2)
signal.seas.2   <- sigex.signal(data,param,mdl,3)
signal.seas.3   <- sigex.signal(data,param,mdl,4)
signal.seas.4   <- sigex.signal(data,param,mdl,5)
signal.seas.5   <- sigex.signal(data,param,mdl,6)
signal.seas     <- sigex.signal(data,param,mdl,2:6)
signal.sa       <- sigex.signal(data,param,mdl,c(1,7))

extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
extract.seas     <- sigex.extract(data,signal.seas,mdl,param)
extract.seas.1   <- sigex.extract(data,signal.seas.1,mdl,param)
extract.seas.2   <- sigex.extract(data,signal.seas.2,mdl,param)
extract.seas.3   <- sigex.extract(data,signal.seas.3,mdl,param)
extract.seas.4   <- sigex.extract(data,signal.seas.4,mdl,param)
extract.seas.5   <- sigex.extract(data,signal.seas.5,mdl,param)
extract.sa       <- sigex.extract(data,signal.sa,mdl,param)

# ---- No band plots ----------------------------------------------------------

subseries <- 1

xss = data[, subseries]
s1  = extract.trendann[[2]][, subseries]
s2  = extract.seas[[2]][, subseries]
s0 = xss - s1 - s2 + 0
s10 = extract.sa[[2]][, subseries]

#modify values for plotting
ms2 = s2 + 30
ms0 = s0 - 20

plot(range(t), range(xss,s1, ms2, ms0), type="n")
lines(as.numeric(xss))
lines(s1, col="tomato")
lines(ms2, col="seagreen")
lines(ms0, col="navyblue")


# ---- FRF --------------------------------------------------------------------

sigex.frf(data = xss, param = param, mdl = mdl, sigcomps = 1, grid = 6000)

# --- Filter weights ----------------------------------------------------------

F = signal.trendann[[1]]
fw = F[250, 1:500]
plot(fw, type="l")
sum(fw)
