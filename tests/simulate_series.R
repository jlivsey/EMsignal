
n=300
Phi=diag(3)
Sig=toeplitz(.9^(0:2))

s1 = gen_trendComp(n, Phi, Sig)
s2 = gen_seasComp(n, Phi, Sig)
s0 = rmvnorm(n = n, mean = rep(0,Ndim), sigma = diag(3))

x = s1+s2+s0

plot(x[,1], type="l")
spec.ar(x[,1])
spec.ar(s2[,1])
spec.ar(s1[,1])

plot(ts(x))


