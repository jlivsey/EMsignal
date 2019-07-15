# source("sim1-B12.R")

d = 12
Gam1 = toeplitz(ARMA2acf(ma=rep(1,11), lag.max = (TT-1-d)))
Gam2 = toeplitz(ARMA2acf(ma=-1, lag.max = (TT-1-d)))
Gam3 = toeplitz(ARMA2acf(ma=c(rep(0,11), -1), lag.max = (TT-1-d)))
Gam  = list(Gam1, Gam2, Gam3)
invGam = lapply(Gam, solve)

# ---- Initialize values for first iteration ----------------------------------

# Choices: par.default, param.mom

# Set Sig
sig2lik(param2sig(par.default))
sig2lik(param2sig(param.mom))
Sig <- param2sig(par.default)

# Set Error matrix
M1 = block2array(signal.trendann[[2]], N = N, TT = TT)
M2 = block2array(signal.seas[[2]],     N = N, TT = TT)
M3 = block2array(signal.irr[[2]],      N = N, TT = TT)
M = list(M1, M2, M3)

# Set signal estimates
S1 = extract.trendann[[1]]
S1d = diff(S1, 12)
S2 = extract.seas[[1]]
S2d = diff(S2, 12)
S3 = extract.irr[[1]]
S3d = diff(S3, 12)
S = list(S1d, S2d, S3d)

lMS = list(M, S)

# -----------------------------------------------------------------------------

iters <- 5
Nc <- length(unlist(Sig)) # number columns of save matrix
Sig.save <- matrix(NA, nrow = iters+1, ncol= Nc+1) # storage container
Sig.save[1, ] <- c(unlist(Sig), sig2lik(Sig)) # first row initial conditions
for(i in 1:iters) {
  out = EMiterate_1_B12(Sig, lMS, data, mdl)
  Sig = out[[1]]
  lMS = out[[2]]
  lik <- sig2lik(Sig)
  print(i)
  print(lik)
  print("------------------------------------")
}


# ---- Plot parameter estimates over time -------------------------------------
library(ggplot2)
library(dplyr)

dat <- data.frame(Sig.save)
colnames(dat) <- c('t11','t21','t12','t22',
                   's11','s21','s12','s22',
                   'i11','i21','i12','i22','lik')
dat$iter <- 1:dim(dat)[1]

# Plot Likelihood
dat.gath <- dat %>%
  select(lik, iter) %>%
  tidyr::gather('variable', 'value', -iter)

ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
  geom_line() +
  geom_hline(yintercept = lik.true)

# Plot all parameters
dat.gath <- dat %>%
  select(-lik) %>%
  tidyr::gather('variable', 'value', -iter)

ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
  geom_line()
