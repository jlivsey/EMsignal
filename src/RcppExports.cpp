// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// subBlock
arma::mat subBlock(const arma::mat& m, const int& N, const int& el, const int& k);
RcppExport SEXP _EMsigex_subBlock(SEXP mSEXP, SEXP NSEXP, SEXP elSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type el(elSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(subBlock(m, N, el, k));
    return rcpp_result_gen;
END_RCPP
}
// matrixDiff
arma::field<arma::mat> matrixDiff(const arma::mat& m, const int& N, const int& TT, const arma::vec& delta);
RcppExport SEXP _EMsigex_matrixDiff(SEXP mSEXP, SEXP NSEXP, SEXP TTSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixDiff(m, N, TT, delta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EMsigex_subBlock", (DL_FUNC) &_EMsigex_subBlock, 4},
    {"_EMsigex_matrixDiff", (DL_FUNC) &_EMsigex_matrixDiff, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_EMsigex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
