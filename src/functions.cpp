// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// [[Rcpp::export]]
arma::mat subBlock(const arma::mat& m, int N, int el, int k){
  // This function assumes R indexing is being passed
  return m.submat(N*(el-1), N*(k-1), N*el-1, N*k-1);
}

// [[Rcpp::export]]
arma::field<arma::mat> matrixDiff(const arma::mat& m, const int& N,
                                  const int& TT, const arma::vec& delta) {
  // convert input block matrix to 4-touple field
  arma::field<arma::mat> M(TT,TT);
  for (int el = 0; el < TT; el++) {
    for (int k = 0; k < TT; k++) {
      M(el, k) = subBlock(m, N, el+1, k+1);
    }
  }
  //
  int dd = delta.size() - 1;
  // Rprintf("dd = %i", dd);
  arma::field<arma::mat> Mout(TT-dd, TT-dd);
  arma::mat Z(N,N);
  Z.fill(0);
  for (int el = 0; el < (TT-dd); el++) {
    for (int k = 0; k < (TT-dd); k++) {
      Mout(el, k) = Z;
      for (int t = 0; t < (dd+1); t++ ){
        for (int s = 0; s < (dd+1); s++ ){
// Rprintf("el=%i, k=%i, t=%i, s=%i \n", el, k, t, s);
          Mout(el, k) += delta(t) * delta(s) * M(el-t+dd, k-s+dd);
        }
      }
    }
  }
  return Mout;
}
