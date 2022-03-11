library(sigex)

data = matrix(c(retail$`Clothing stores.dat`,
                retail$`Jewelry stores.dat`),
              ncol=2, nrow = 282, byrow = FALSE)
data = demean(data)


N = 2
T <- TT <- 282
t = 1:T


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

# ---- MLE --------------------------------------------------------------------
# MLE estimation of full model

# setup, and fix parameters as desired
mdl.mle <- mdl
psi.mle <- psi.default
flag.mle <- Im(psi.mle)
par.mle <- sigex.psi2par(psi.mle,mdl.mle,data)

# run fitting
fit.mle <- sigex.mlefit(data,par.mle,flag.mle,mdl.mle,"bfgs")

psi.mle[flag.mle==1] <- fit.mle[[1]]$par
psi.mle <- psi.mle + 1i*flag.mle
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]
resid.mle <- sigex.resid(psi.mle,mdl.mle,data)
acf(t(resid.mle),lag.max=40)

print(eigen(hess)$values)
taus <- log(sigex.conditions(data,psi.mle,mdl.mle))
print(taus)

tstats <- sigex.tstats(mdl.mle,psi.mle,hess)
stderrs <- sigex.psi2par(tstats,mdl,data)
print(tstats)

# bundle for default span
analysis.mle <- sigex.bundle(data,transform,mdl.mle,psi.mle)

psi <- analysis.mle[[4]]
param.mle <- sigex.psi2par(psi,mdl,data)


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
fw = FF[1, 2, T/2, ]
plot(fw, type="l")
abline(v=seq(T/2, T, 12), lty="dotted")
sum(fw)
