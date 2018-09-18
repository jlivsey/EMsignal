library(sigex)
library(mvtnorm)

# ---- Simulate Data ----------------------------------------------------------
N = 3
T <- TT <- 300
t = 1:T
Phi=diag(N)
Sig=diag(N)

s1 = gen_trendComp(T, Phi, Sig)
s2 = gen_seasComp(T, Phi, Sig)
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
data = demean(data)
plot(ts(data))

spec.ar(data[,1])
abline(v=seq(1/12, 1/2, 1/12))

x12 = diff(data[,1], 12)
mean(x12)

# ---- Modeling ---------------------------------------------------------------
transform = "none"
x <- t(data)
N <- dim(x)[1]
T <- dim(x)[2]

# ---- Model ------------------------------------------------------------------
# stochastic effects
def <- c(0,1,0,1)

mdl <- NULL
# 1 trend-cycle (annual cycle)
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-1),def)
# 2 first atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(2*pi/12),1),def)
# 3 second atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(4*pi/12),1),def)
# 4 third atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(6*pi/12),1),def)
# 5 fourth atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(8*pi/12),1),def)
# 6 fifth atomic monthly seas
mdl <- sigex.add(mdl,seq(1,N),"wn",c(1,-2*cos(10*pi/12),1),def)
# 7 pi freq
mdl <- sigex.add(mdl,seq(1,N),"wn", c(1,1), def)
# 8 irregular
mdl <- sigex.add(mdl,seq(1,N),"wn",1,def)
# fixed effects [mean effect]
mdl <- sigex.meaninit(mdl,data,0)

# ---- Model 2 ----------------------------------------------------------------
mdl = NULL
mdl = sigex.add(mdl, seq(1,N), 'wn', c(1,-1), def)
mdl = sigex.add(mdl, seq(1,N), 'wn', rep(1, 12), def)
mdl <- sigex.add(mdl,seq(1,N),"wn",1,def)
mdl <- sigex.meaninit(mdl,data,0)


# Set default parameters
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


# ---- DMA for model 1 --------------------------------------------------------
signal.trendann <- sigex.signal(data,param,mdl,1)
signal.seas.1   <- sigex.signal(data,param,mdl,2)
signal.seas.2   <- sigex.signal(data,param,mdl,3)
signal.seas.3   <- sigex.signal(data,param,mdl,4)
signal.seas.4   <- sigex.signal(data,param,mdl,5)
signal.seas.5   <- sigex.signal(data,param,mdl,6)
signal.seas.6   <- sigex.signal(data,param,mdl,7)
signal.seas     <- sigex.signal(data,param,mdl,2:7)
signal.sa       <- sigex.signal(data,param,mdl,c(1,8))
signal.irr      <- sigex.signal(data,param,mdl,8)

extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
extract.seas     <- sigex.extract(data,signal.seas,mdl,param)
extract.seas.1   <- sigex.extract(data,signal.seas.1,mdl,param)
extract.seas.2   <- sigex.extract(data,signal.seas.2,mdl,param)
extract.seas.3   <- sigex.extract(data,signal.seas.3,mdl,param)
extract.seas.4   <- sigex.extract(data,signal.seas.4,mdl,param)
extract.seas.5   <- sigex.extract(data,signal.seas.5,mdl,param)
extract.sa       <- sigex.extract(data,signal.sa,mdl,param)
extract.irr      <- sigex.extract(data,signal.irr,mdl,param)

# ---- DMA for model 2 --------------------------------------------------------
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
fw = FF[1, 1, T/2, ]
plot(fw, type="l")
abline(v=seq(T/2, T, 12), lty="dotted")
sum(fw)

par(mfrow=c(3,2), mar=c(2,2,3,1))
for(i in 1:N){
  for(j in i:N){
    cat(i,j)
    fw = FF[i, j, T/2, ]
    plot(fw, type="l", main=paste(i,j))
    abline(v=seq(T/2, T, 12), lty="dotted")
  }
}
