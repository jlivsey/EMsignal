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

M1 = block2array(signal.trendann[[2]], N = N, TT = TT) # M1_t = M1[, t, , t]
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

iters <- 30
Nc <- length(unlist(Sig.mom))
Sig.save <- matrix(NA, nrow = iters+2, ncol= Nc+1)
Sig.save[1, ] <- c(unlist(Sig.true), lik.true)
Sig.save[2, ] <- c(unlist(Sig.mom), lik.mom)
for(i in 16:iters) {
  if(i==1) Sig <- Sig.mom
  out = EMiterate_1_B12(Sig, lMS, data, mdl, invGam)
  Sig = out[[1]]
  lMS = out[[2]]
  lik <- sig2lik(Sig, mdl, data)
  Sig.save[i+2, ] <- c(unlist(Sig), lik)

  print(Sys.time())
  print(i)
  print(lik)
  print("------------------------------------")
}

# MOM estimates
Sig_trend_mom = Sig.save[2, 1:9]   |> matrix(3, 3)
Sig_seas_mom  = Sig.save[2, 10:18] |> matrix(3, 3)
Sig_irr_mom   = Sig.save[2, 19:27] |> matrix(3, 3)

# EM estimates
Sig_trend_em = Sig.save[iters, 1:9]   |> matrix(3, 3)
Sig_seas_em  = Sig.save[iters, 10:18] |> matrix(3, 3)
Sig_irr_em   = Sig.save[iters, 19:27] |> matrix(3, 3)

# Compare
cbind(Sig_trend_mom, Sig_trend_em, Sig_trend_mom - Sig_trend_em) |> round(3)
cbind(Sig_seas_mom, Sig_seas_em, Sig_seas_mom - Sig_seas_em) |> round(3)
cbind(Sig_irr_mom, Sig_irr_em, Sig_irr_mom - Sig_irr_em) |> round(3)



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













