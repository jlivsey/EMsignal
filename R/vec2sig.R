

vec2sig <- function(vec, N){
  numSigs <- length(vec)/N^2
  Sig <- vector(mode = 'list', length = numSigs)
  for(i in 1:numSigs){
    dat <- vec[(N^2*(i-1)+1):(N^2*i)]
    Sig[[i]] <- matrix(dat, N, N)
  }
  return(Sig)
}
