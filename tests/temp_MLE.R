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

