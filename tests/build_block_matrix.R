
N = 3
TT = 300

Idx = matrix(1:(N*TT), ncol=N, byrow = TRUE)
# Idx

A = matrix(NA, N*TT, N*TT)
counter = 0
for(i in 1:dim(Idx)[1]){
  for(j in 1:dim(Idx)[1]){
    A[Idx[i, ], Idx[j,]] = counter
    counter = counter + 1
  }
}
