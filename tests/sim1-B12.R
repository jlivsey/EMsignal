library(sigex)
library(mvtnorm)

# ---- Simulate Data ----------------------------------------------------------
N = 2
T <- TT <- 144
t = 1:T
Phi=diag(N)
Sig=diag(N)

s1 = gen_trendComp(T, Phi, Sig)
s2 = gen_seasComp(T, Phi, Sig)
s0 = rmvnorm(n = T, mean = rep(0,N), sigma = diag(N))

data = s1+s2+s0
data = demean(data)
plot(ts(data))

# ---- Modeling ---------------------------------------------------------------
transform = "none"
x <- t(data)
N <- dim(x)[1]
T <- dim(x)[2]

# ---- Model ------------------------------------------------------------------
def <- c(0,1,0,1)
mdl = NULL
mdl = sigex.add(mdl, seq(1,N), 'wn', c(1,-1), def)    # Trend
mdl = sigex.add(mdl, seq(1,N), 'wn', rep(1, 12), def) # Seasonal
mdl <- sigex.add(mdl,seq(1,N),"wn",1,def)             # Irregular
mdl <- sigex.meaninit(mdl,data,0)

# Set default parameters
par.default <- sigex.default(mdl,data)[[1]]
flag.default <- sigex.default(mdl,data)[[2]]
psi.default <- sigex.par2psi(par.default,flag.default,mdl)

# Set param to TRUE values
param = par.default

# ---- DMA for model ----------------------------------------------------------
signal.trendann <- sigex.signal(data,param,mdl,1)
signal.seas     <- sigex.signal(data,param,mdl,2)
signal.irr      <- sigex.signal(data,param,mdl,3)

extract.trendann <- sigex.extract(data,signal.trendann,mdl,param)
extract.seas     <- sigex.extract(data,signal.seas,mdl,param)
extract.irr      <- sigex.extract(data,signal.irr,mdl,param)


# ---- Plots ------------------------------------------------------------------
subseries <- 1
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
abline(v=seq(0, T/2, 12), lty="dotted")
sum(fw)


