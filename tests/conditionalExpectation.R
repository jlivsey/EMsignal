

conditionalExpectiation <- function(){
  out = 0
  for(j in J){

    s = sigex.extract[[]][]
    M =

    out = out + tr(M + s %*% t(s))
  }
  return(out)
}


for(j in 1:J){

  lsigma = list()

  val = 0
  for(k in 1:T){
    new.val = lsignal[[j]][[2]][k, k] + lextract[[j]][[1]][k, ]
    val = val + new.val
  }
  val = val / T

}


