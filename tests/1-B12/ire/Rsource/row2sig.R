#' convert vector of stored values to Sig list form
#'
#' @param v vector of unlisted values
#' @param N dimension of output matricies
#'
#' @return list of Sig matricies
#' @export
#'

row2sig <- function(v, N){
  Sig <- list()
  J <- floor(length(v)/N^2)
  for(i in 1:J){
    idx <- (i-1) * N^2 + 1
    Sig[[i]] <- matrix(v[idx:(idx + (N^2-1))], N, N)
  }
  return(Sig)
}
