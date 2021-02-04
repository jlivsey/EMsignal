d = 12
Gam1 = toeplitz(ARMAauto(ma = rep(1,11), ar = NULL, lag.max = (TT-1-d)))
Gam2 = toeplitz(ARMAauto(ma = -1, ar = NULL, lag.max = (TT-1-d)))
Gam3 = toeplitz(ARMAauto(ma = c(rep(0,11), -1), ar = NULL, lag.max = (TT-1-d)))
Gam  = list(Gam1, Gam2, Gam3)
invGam = lapply(Gam, solve)

# ---- Initialize values for first iteration ----------------------------------

# Choices: par.default, param.mom

# Set Sig
sig2lik(param2sig(par.default), mdl, data)
sig2lik(param2sig(param.mom),   mdl, data)
Sig <- param2sig(param.mom)

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

# ---- Main EM loop ------------------------------------------------------------

iters <- 20
Nc <- length(unlist(Sig)) # number columns of save matrix
Sig.save <- matrix(NA, nrow = iters+1, ncol= Nc+1) # storage container
Sig.save[1, ] <- c(unlist(Sig), sig2lik(Sig, mdl, data)) # first row initial conditions
for(i in 11:iters) {
  out = EMiterate_1_B12(Sig, lMS, data, mdl, invGam)
  Sig = out[[1]]
  lMS = out[[2]]
  lik <- sig2lik(Sig, mdl , data)
  print(i)
  print(lik)
  print("------------------------------------")
  Sig.save[i+1, ] <- c(unlist(Sig), lik)
  save(Sig.save, file = "20210203-Sigsave.Rdata")
}


# ---- Plot parameter estimates over time -------------------------------------
library(ggplot2)
library(dplyr)

dat <- data.frame(Sig.save)
dat$iter <- 1:dim(dat)[1]
dat2 <- dat %>%
  rename(lik = X49)

# Plot Likelihood
dat.gath <- dat2 %>%
  tidyr::gather('variable', 'value', -iter, -lik)

ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
  geom_line()
  # geom_hline(yintercept = lik.true)

# Plot all parameters
dat.gath <- dat %>%
  select(-lik) %>%
  tidyr::gather('variable', 'value', -iter)

ggplot(dat.gath, aes(x = iter, y = value, colour = variable)) +
  geom_line()
