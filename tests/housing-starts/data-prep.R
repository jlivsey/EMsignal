library(sigex)
# library(mvtnorm)
library(EMsigex)

# ---- Prep ---------------------------------------------------------------
start.date = c(1964,1)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(starts,start.date,period,c("South","West","NE","MW"),TRUE)

#############################
## select span and transforms

## all data for NE-MW with log transform
transform <- "log"
aggregate <- FALSE
subseries <- c(1,2,3,4)
begin.date <- start(dataALL.ts)
end.date <- end(dataALL.ts)
range <- NULL
data <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)

N <- dim(starts)[2]
T <- dim(starts)[1]
TT <- T

# ---- Model ------------------------------------------------------------------
mdl <- NULL
# mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-1))   # first diff
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"trend",c(1,-2,1)) # second diff
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"seasonal", rep(1,12))
mdl <- sigex.add(mdl,seq(1,N),"arma",c(0,0),0,"irregular",1)
# regressors:
mdl <- sigex.meaninit(mdl, data, 0)

# Set default parameters
constraint <- NULL
par.default <- sigex.default(mdl, data, constraint)
psi.default <- sigex.par2psi(par.default, mdl)

# Set param to TRUE values
param = par.default

# # ---- MOM estimates for param ------------------------------------------------
mdl.mom <- mdl
par.mom <- sigex.momfit(data, par.default, mdl.mom)
psi.mom <- sigex.par2psi(par.mom, mdl.mom)
#resid.mom <- sigex.resid(psi.mom,mdl.mom,data)

thresh <- -1.66

if(N > 1) {
  reduced.mom <- sigex.reduce(data, par.mom, mdl.mom, thresh, FALSE)
  mdl.mom <- reduced.mom[[1]]
  par.mom <- reduced.mom[[2]]
  psi.mom <- sigex.par2psi(par.mom, mdl.mom)
  #resid.mom <- sigex.resid(psi.mom,mdl.mom,data)
}

# bundle for default span
analysis.mom <- sigex.bundle(data, transform, mdl.mom, psi.mom)

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

# pdf(file = "south-EM.pdf", width = 6, height = 4)
subseries <- 1
xss = data[, subseries]
s1.hat  = extract.trendann[[1]][, subseries]
s2.hat  = extract.seas[[1]][, subseries]
s0.hat = extract.irr[[1]][, subseries]
{
  op = par(mfrow=c(3,1), mar=c(2,3,2,1))
  plot(as.numeric(xss), type="l", main = "South")
  lines(s1.hat, col="tomato")
  plot(s2.hat, type="l", col="seagreen"); abline(h=0, lty="dotted")
  abline(v=seq(1,TT,12), lty="dashed")
  plot(s0.hat, type="l", col="navyblue"); abline(h=0, lty="dotted")
  par(op)
}
# dev.off()

# --- Filter weights ----------------------------------------------------------
FF = signal.trendann[[1]]
FF = block2array(FF, N, T)
fw = FF[1, T/2, 1, ]
plot(fw, type="l")
abline(v=seq(T/2, T, 12), lty="dotted")
sum(fw)

