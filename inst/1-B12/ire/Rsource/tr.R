#' Calculate trace of a square matrix
#'
#' @param m square matrix to calculate trace of
#'
#' @return trace of given matrix
#' @export
#'

tr <- function (m)
{
  total_sum <- 0
  if(is.matrix(m))
  {
    row_count <- nrow(m)
    col_count <- ncol(m)
    if(row_count == col_count)
    {
      total_sum <-sum(diag(m))
      total_sum
    }
    else
    {
      message ('Matrix is not square')
    }
  }
  else
  {
    message( 'Object is not a matrix')

  }
}
