// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// DST_Solver
Eigen::VectorXd DST_Solver(int Dim, const Eigen::VectorXd& Lambda, const Eigen::SparseMatrix<double>& B, const Eigen::VectorXd& U, int Ns);
RcppExport SEXP _YauYauAL_DST_Solver(SEXP DimSEXP, SEXP LambdaSEXP, SEXP BSEXP, SEXP USEXP, SEXP NsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Dim(DimSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    rcpp_result_gen = Rcpp::wrap(DST_Solver(Dim, Lambda, B, U, Ns));
    return rcpp_result_gen;
END_RCPP
}
// ExpandGrid
Rcpp::NumericMatrix ExpandGrid(int Dim, Rcpp::NumericVector s);
RcppExport SEXP _YauYauAL_ExpandGrid(SEXP DimSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Dim(DimSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(ExpandGrid(Dim, s));
    return rcpp_result_gen;
END_RCPP
}
// NormalizedExp
arma::vec NormalizedExp(const arma::vec& x);
RcppExport SEXP _YauYauAL_NormalizedExp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(NormalizedExp(x));
    return rcpp_result_gen;
END_RCPP
}
// computeB
Eigen::SparseMatrix<double> computeB(const Eigen::MatrixXd& s, const std::vector<Eigen::SparseMatrix<double>>& D_list, double Dt, double Ds, Rcpp::Function f, Rcpp::Function df, Rcpp::Function h);
RcppExport SEXP _YauYauAL_computeB(SEXP sSEXP, SEXP D_listSEXP, SEXP DtSEXP, SEXP DsSEXP, SEXP fSEXP, SEXP dfSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::SparseMatrix<double>>& >::type D_list(D_listSEXP);
    Rcpp::traits::input_parameter< double >::type Dt(DtSEXP);
    Rcpp::traits::input_parameter< double >::type Ds(DsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type df(dfSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(computeB(s, D_list, Dt, Ds, f, df, h));
    return rcpp_result_gen;
END_RCPP
}
// computeLambda
Eigen::VectorXd computeLambda(int Dim, int Ns, double Dt, double Ds);
RcppExport SEXP _YauYauAL_computeLambda(SEXP DimSEXP, SEXP NsSEXP, SEXP DtSEXP, SEXP DsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Dim(DimSEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< double >::type Dt(DtSEXP);
    Rcpp::traits::input_parameter< double >::type Ds(DsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLambda(Dim, Ns, Dt, Ds));
    return rcpp_result_gen;
END_RCPP
}
// generateD
std::vector<Eigen::SparseMatrix<double>> generateD(int Dim, int Ns, double Ds);
RcppExport SEXP _YauYauAL_generateD(SEXP DimSEXP, SEXP NsSEXP, SEXP DsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Dim(DimSEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< double >::type Ds(DsSEXP);
    rcpp_result_gen = Rcpp::wrap(generateD(Dim, Ns, Ds));
    return rcpp_result_gen;
END_RCPP
}
// Simulate_State_Obser
Rcpp::List Simulate_State_Obser(double Dt, int Ntau, int NtNtau, Rcpp::Function f, Rcpp::Function h, int Dim, Rcpp::Nullable<int> seed);
RcppExport SEXP _YauYauAL_Simulate_State_Obser(SEXP DtSEXP, SEXP NtauSEXP, SEXP NtNtauSEXP, SEXP fSEXP, SEXP hSEXP, SEXP DimSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Dt(DtSEXP);
    Rcpp::traits::input_parameter< int >::type Ntau(NtauSEXP);
    Rcpp::traits::input_parameter< int >::type NtNtau(NtNtauSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type Dim(DimSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_YauYauAL_DST_Solver", (DL_FUNC) &_YauYauAL_DST_Solver, 5},
    {"_YauYauAL_ExpandGrid", (DL_FUNC) &_YauYauAL_ExpandGrid, 2},
    {"_YauYauAL_NormalizedExp", (DL_FUNC) &_YauYauAL_NormalizedExp, 1},
    {"_YauYauAL_computeB", (DL_FUNC) &_YauYauAL_computeB, 7},
    {"_YauYauAL_computeLambda", (DL_FUNC) &_YauYauAL_computeLambda, 4},
    {"_YauYauAL_generateD", (DL_FUNC) &_YauYauAL_generateD, 3},
    {"_YauYauAL_Simulate_State_Obser", (DL_FUNC) &_YauYauAL_Simulate_State_Obser, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_YauYauAL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
