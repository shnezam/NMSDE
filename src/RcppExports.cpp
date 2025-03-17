// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// permute
arma::cube permute(arma::cube a, int order);
RcppExport SEXP _NMSDE_permute(SEXP aSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(permute(a, order));
    return rcpp_result_gen;
END_RCPP
}
// inChol2spec_c
arma::mat inChol2spec_c(arma::mat inchol, List init);
RcppExport SEXP _NMSDE_inChol2spec_c(SEXP incholSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type inchol(incholSEXP);
    Rcpp::traits::input_parameter< List >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(inChol2spec_c(inchol, init));
    return rcpp_result_gen;
END_RCPP
}
// tA2G_c_cube
arma::cube tA2G_c_cube(arma::cube thetas, arma::cube As, List init);
RcppExport SEXP _NMSDE_tA2G_c_cube(SEXP thetasSEXP, SEXP AsSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type As(AsSEXP);
    Rcpp::traits::input_parameter< List >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(tA2G_c_cube(thetas, As, init));
    return rcpp_result_gen;
END_RCPP
}
// loglike_c_cube
Rcpp::List loglike_c_cube(arma::cube thetas, arma::cube As, arma::vec lambda, Rcpp::List init, Rcpp::List tildeX);
RcppExport SEXP _NMSDE_loglike_c_cube(SEXP thetasSEXP, SEXP AsSEXP, SEXP lambdaSEXP, SEXP initSEXP, SEXP tildeXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type As(AsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type init(initSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tildeX(tildeXSEXP);
    rcpp_result_gen = Rcpp::wrap(loglike_c_cube(thetas, As, lambda, init, tildeX));
    return rcpp_result_gen;
END_RCPP
}
// Workingfish_cube
Rcpp::List Workingfish_cube(arma::cube thetas, arma::cube As, arma::vec lambda, Rcpp::List init, Rcpp::List tildeX);
RcppExport SEXP _NMSDE_Workingfish_cube(SEXP thetasSEXP, SEXP AsSEXP, SEXP lambdaSEXP, SEXP initSEXP, SEXP tildeXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type As(AsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type init(initSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tildeX(tildeXSEXP);
    rcpp_result_gen = Rcpp::wrap(Workingfish_cube(thetas, As, lambda, init, tildeX));
    return rcpp_result_gen;
END_RCPP
}
// tA2G_c_mat
arma::cube tA2G_c_mat(arma::mat thetas, arma::cube As, Rcpp::List init);
RcppExport SEXP _NMSDE_tA2G_c_mat(SEXP thetasSEXP, SEXP AsSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type As(AsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(tA2G_c_mat(thetas, As, init));
    return rcpp_result_gen;
END_RCPP
}
// loglike_c_mat
Rcpp::List loglike_c_mat(arma::mat thetas, arma::cube As, arma::vec lambda, List init, Rcpp::List tildeX);
RcppExport SEXP _NMSDE_loglike_c_mat(SEXP thetasSEXP, SEXP AsSEXP, SEXP lambdaSEXP, SEXP initSEXP, SEXP tildeXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type As(AsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type init(initSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tildeX(tildeXSEXP);
    rcpp_result_gen = Rcpp::wrap(loglike_c_mat(thetas, As, lambda, init, tildeX));
    return rcpp_result_gen;
END_RCPP
}
// Workingfish_mat
Rcpp::List Workingfish_mat(arma::mat thetas, arma::cube As, arma::vec lambda, Rcpp::List init, Rcpp::List tildeX);
RcppExport SEXP _NMSDE_Workingfish_mat(SEXP thetasSEXP, SEXP AsSEXP, SEXP lambdaSEXP, SEXP initSEXP, SEXP tildeXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type As(AsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type init(initSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tildeX(tildeXSEXP);
    rcpp_result_gen = Rcpp::wrap(Workingfish_mat(thetas, As, lambda, init, tildeX));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NMSDE_permute", (DL_FUNC) &_NMSDE_permute, 2},
    {"_NMSDE_inChol2spec_c", (DL_FUNC) &_NMSDE_inChol2spec_c, 2},
    {"_NMSDE_tA2G_c_cube", (DL_FUNC) &_NMSDE_tA2G_c_cube, 3},
    {"_NMSDE_loglike_c_cube", (DL_FUNC) &_NMSDE_loglike_c_cube, 5},
    {"_NMSDE_Workingfish_cube", (DL_FUNC) &_NMSDE_Workingfish_cube, 5},
    {"_NMSDE_tA2G_c_mat", (DL_FUNC) &_NMSDE_tA2G_c_mat, 3},
    {"_NMSDE_loglike_c_mat", (DL_FUNC) &_NMSDE_loglike_c_mat, 5},
    {"_NMSDE_Workingfish_mat", (DL_FUNC) &_NMSDE_Workingfish_mat, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_NMSDE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
