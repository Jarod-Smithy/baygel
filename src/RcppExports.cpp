// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// blockBAGENI
List blockBAGENI(arma::mat X, int burnin, int iterations, bool verbose, double r, double s, double a, double b);
RcppExport SEXP _baygel_blockBAGENI(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP verboseSEXP, SEXP rSEXP, SEXP sSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBAGENI(X, burnin, iterations, verbose, r, s, a, b));
    return rcpp_result_gen;
END_RCPP
}
// blockBAGENII
List blockBAGENII(arma::mat X, int burnin, int iterations, bool verbose, double s, double b);
RcppExport SEXP _baygel_blockBAGENII(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP verboseSEXP, SEXP sSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBAGENII(X, burnin, iterations, verbose, s, b));
    return rcpp_result_gen;
END_RCPP
}
// blockBAGL
List blockBAGL(arma::mat X, int burnin, int iterations, bool verbose, double r, double s);
RcppExport SEXP _baygel_blockBAGL(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP verboseSEXP, SEXP rSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBAGL(X, burnin, iterations, verbose, r, s));
    return rcpp_result_gen;
END_RCPP
}
// blockBAGR
List blockBAGR(arma::mat X, int burnin, int iterations, bool verbose, double a, double b);
RcppExport SEXP _baygel_blockBAGR(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP verboseSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBAGR(X, burnin, iterations, verbose, a, b));
    return rcpp_result_gen;
END_RCPP
}
// blockBGEN
List blockBGEN(arma::mat X, int burnin, int iterations, double lambda, double sig, bool verbose);
RcppExport SEXP _baygel_blockBGEN(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP lambdaSEXP, SEXP sigSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBGEN(X, burnin, iterations, lambda, sig, verbose));
    return rcpp_result_gen;
END_RCPP
}
// blockBGL
List blockBGL(arma::mat X, int burnin, int iterations, double lambda, bool verbose);
RcppExport SEXP _baygel_blockBGL(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP lambdaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBGL(X, burnin, iterations, lambda, verbose));
    return rcpp_result_gen;
END_RCPP
}
// blockBGR
List blockBGR(arma::mat X, int burnin, int iterations, double sig, bool verbose);
RcppExport SEXP _baygel_blockBGR(SEXP XSEXP, SEXP burninSEXP, SEXP iterationsSEXP, SEXP sigSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(blockBGR(X, burnin, iterations, sig, verbose));
    return rcpp_result_gen;
END_RCPP
}
