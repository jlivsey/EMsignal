d.full = length(diff.full)
# Build differencing matrix
delta0pad = c(diff.full, rep(0, TT-d.full+1))
D = suppressWarnings(matrix(delta0pad, nrow = TT-d.full,
                            ncol = TT, byrow = TRUE))

delta = c(1, -1)
TT = 10

delta0pad = c(delta, rep(0, TT-length(delta)+1))

D = suppressWarnings(matrix(delta0pad, TT-length(delta)+1, TT, byrow=TRUE))

x = cbind(1:10, 1:10)
x
D%*%x
