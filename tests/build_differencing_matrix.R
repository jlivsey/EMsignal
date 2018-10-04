d.full = length(diff.full)
# Build differencing matrix
delta0pad = c(diff.full, rep(0, TT-d.full+1))
D = suppressWarnings(matrix(delta0pad, nrow = TT-d.full,
                            ncol = TT, byrow = TRUE))
