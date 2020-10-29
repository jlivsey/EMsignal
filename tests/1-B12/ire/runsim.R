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

# Load
start.date <- c(1990, 1)
period <- 12
dataALL.ts <- sigex.load(data = data,
                         start.date = start.date,
                         period = period,
                         epithets = c("Dim1","Dim2","Dim3"),
                         plot = TRUE)

# Prep for sigex
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2,3)
begin.date <- start.date
end.date <- end(dataALL.ts)
range <- list(begin.date,end.date)
data.ts <- sigex.prep(data.ts = dataALL.ts,
                      transform = transform,
                      aggregate = aggregate,
                      subseries = subseries,
                      range = range,
                      plot = TRUE)

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
mdl <- sigex.meaninit(mdl, data.ts, 0)

# Set default parameters
constraint <- NULL
par.default <- sigex.default(mdl, data.ts, constraint)
psi.default <- sigex.par2psi(par.default, mdl)

# Set param to TRUE values
param = par.default

# # ---- MOM estimates for param ------------------------------------------------
mdl.mom <- mdl
par.mom <- sigex.momfit(data,par.default,mdl.mom)
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
analysis.mom <- sigex.bundle(data.ts ,transform,mdl.mom,psi.mom)
psi <- analysis.mom[[4]]
param.mom <- sigex.psi2par(psi,mdl,data.ts)

# ---- MLE --------------------------------------------------------------------
#############
# MLE Fitting

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
                        debug  = TRUE)

# set psi from mle estimation
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian

# bundle for default span
analysis.mle <- sigex.bundle(data.ts, transform, mdl, psi.mle)
par.mle <- sigex.psi2par(psi.mle, mdl, data.ts)

# ---- signal extraction from fit - Direct Matrix Approach ---------------------

# which params to use default/mom/mle
param <- par.mle

signal.trendann <- sigex.signal(data.ts, param, mdl, 1)
signal.seas     <- sigex.signal(data.ts, param, mdl, 2)
signal.irr      <- sigex.signal(data.ts, param, mdl, 3)

extract.trendann <- sigex.extract(data.ts, signal.trendann, mdl, param)
extract.seas     <- sigex.extract(data.ts, signal.seas,     mdl, param)
extract.irr      <- sigex.extract(data.ts, signal.irr,      mdl, param)


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

source("tests/1-B12/sim_1-B12.R")

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
(lik.true <- sig2lik(Sig.true, mdl, data))


# ---- Initialize values for first iteration ----------------------------------

param = param.mom
Sig.mom = param2sig(param)
(lik.mom <- sig2lik(Sig.mom, mdl, data))

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

# Sig.mle = param2sig(param.mle)

iters <- 5
Nc <- length(unlist(Sig.mom))
Sig.save <- matrix(NA, nrow = iters+2, ncol= Nc+1)
Sig.save[1, ] <- c(unlist(Sig.true), lik.true)
Sig.save[2, ] <- c(unlist(Sig.mom), lik.mom)
for(i in 1:iters) {
  if(i==1) Sig <- Sig.mom
  out = EMiterate_1_B12(Sig, lMS, data, mdl)
  Sig = out[[1]]
  lMS = out[[2]]
  lik <- sig2lik(Sig, mdl, data)
  Sig.save[i+2, ] <- c(unlist(Sig), lik)

  print(Sys.time())
  print(i)
  print(lik)
  print("------------------------------------")
}



# # ---- Plot parameter estimates over time -------------------------------------
library(ggplot2)
library(dplyr)
library(gridExtra)

dat <- data.frame(Sig.save)
dat$iter <- 1:dim(dat)[1]
colnames(dat) <- c('t11','t21', 't31','t12', 't22', 't32', 't13', 't23', 't33',
                   's11','s21', 's31','s12', 's22', 's32', 's13', 's23', 's33',
                   'i11','i21', 'i31','i12', 'i22', 'i32', 'i13', 'i23', 'i33',
                   'lik', 'iter')

# Plot Likelihood
dat.gath <- dat %>%
  select(lik, iter) %>%
  tidyr::gather('variable', 'value', -iter)

gg_lik <-
  ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
  geom_line() +
  geom_hline(yintercept = lik.true)

# Plot all parameters
dat.gath <- dat %>%
  select(-lik) %>%
  tidyr::gather('variable', 'value', -iter)

gg_param <-
  ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
  geom_line()

gridExtra::grid.arrange(gg_param, gg_lik, ncol = 2)


# pdf("EM-sim-results-2020-03-19.pdf")
# for(i in 1:28){
#   plot(ss[, i], type = 'b')
#   abline(h = ss[1, i])
# }
# dev.off()













