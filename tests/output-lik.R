

# Define Sigma for which you want to see likelihood.
Sig4 <- out[[1]]

# Put in param form
param <- sig2param(Sig4)

# Put in psi form
psi <- sigex.par2psi(param, flag.mom, mdl)

# print likelihood
sigex.lik(psi, mdl, data)
