#' Convert block matrix to array
#'
#' Convert NTxNT matrix with NxN sub-blocks into [N,T,N,T] array
#'
#' @param M NTxNT matrix with NxN sub-blocks
#' @param N dim of sub-blocks
#' @param TT Number of sub-blocks in each row (column)
#'
#' @return an array of dim [N,T,N,T] where [, i, , j] gives the (i,j)th NxN sub-block
#' @export
#'

block2array = function(M, N, TT){
  # if(!all(dim(M)==c(N*TT, N*TT))) stop("M is not (NT)x(NT) matrix")
  # Idx = matrix(1:(N*TT), ncol=N, byrow = TRUE)
  # A = array(dim = c(N,N,TT,TT))
  # for(i in 1:TT){
  #   for(j in 1:TT){
  #     A[,,i,j] = M[Idx[i, ], Idx[j, ]]
  #   }
  # }
  # return(A)
  A = array(M, c(N, TT, N, TT))
  return(A)
}

