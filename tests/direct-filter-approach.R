d = c(-1, 1)
TT = 100

# Construct differencing matrix
d.vec = c(d, rep(0,TT-1))
D = matrix(d.vec, nrow = TT-1, ncol = TT, byrow = TRUE)

# Filter matrix F (from Tucker paper)
F = solve(t(D) %*% D + diag(TT))

# plot concurrent and middle of span filter weights
plot(range(1:TT), range(D), type="n")
lines(F[5,], col="red")
lines(F[1,], col="blue")
lines(F[3,], col = "green")
