

conditionalExpectiation <- function(){
  out = 0
  for(j in J){

    s = sigex.extract[[]][]
    M =

    out = out + tr(M + s %*% t(s))
  }
  return(out)
}
