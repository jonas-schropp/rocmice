// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// tpr
double tpr(IntegerVector r, IntegerVector p);
RcppExport SEXP _rocmice_tpr(SEXP rSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(tpr(r, p));
    return rcpp_result_gen;
END_RCPP
}
// fpr
double fpr(IntegerVector r, IntegerVector p);
RcppExport SEXP _rocmice_fpr(SEXP rSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(fpr(r, p));
    return rcpp_result_gen;
END_RCPP
}
// var_prop
NumericVector var_prop(NumericVector est, double n);
RcppExport SEXP _rocmice_var_prop(SEXP estSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type est(estSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(var_prop(est, n));
    return rcpp_result_gen;
END_RCPP
}
// cut_off
IntegerVector cut_off(NumericVector score, double cutoff);
RcppExport SEXP _rocmice_cut_off(SEXP scoreSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(cut_off(score, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// varlogit
NumericVector varlogit(NumericVector var, NumericVector est);
RcppExport SEXP _rocmice_varlogit(SEXP varSEXP, SEXP estSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type var(varSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type est(estSEXP);
    rcpp_result_gen = Rcpp::wrap(varlogit(var, est));
    return rcpp_result_gen;
END_RCPP
}
// apply_cut_off
NumericMatrix apply_cut_off(NumericVector score, NumericVector unique_vals);
RcppExport SEXP _rocmice_apply_cut_off(SEXP scoreSEXP, SEXP unique_valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type unique_vals(unique_valsSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_cut_off(score, unique_vals));
    return rcpp_result_gen;
END_RCPP
}
// count_vals
int count_vals(IntegerVector x, int val);
RcppExport SEXP _rocmice_count_vals(SEXP xSEXP, SEXP valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type val(valSEXP);
    rcpp_result_gen = Rcpp::wrap(count_vals(x, val));
    return rcpp_result_gen;
END_RCPP
}
// apply_metrics
List apply_metrics(const IntegerMatrix& m, const IntegerVector& r);
RcppExport SEXP _rocmice_apply_metrics(SEXP mSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(apply_metrics(m, r));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_rowMeans
NumericVector rcpp_rowMeans(NumericMatrix x);
RcppExport SEXP _rocmice_rcpp_rowMeans(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_rowMeans(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_var
double rcpp_var(NumericVector samples);
RcppExport SEXP _rocmice_rcpp_var(SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_var(samples));
    return rcpp_result_gen;
END_RCPP
}
// rowVars
NumericVector rowVars(NumericMatrix x);
RcppExport SEXP _rocmice_rowVars(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowVars(x));
    return rcpp_result_gen;
END_RCPP
}
// pool_rr_m
List pool_rr_m(NumericMatrix est, NumericMatrix var);
RcppExport SEXP _rocmice_pool_rr_m(SEXP estSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type est(estSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(pool_rr_m(est, var));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rocmice_tpr", (DL_FUNC) &_rocmice_tpr, 2},
    {"_rocmice_fpr", (DL_FUNC) &_rocmice_fpr, 2},
    {"_rocmice_var_prop", (DL_FUNC) &_rocmice_var_prop, 2},
    {"_rocmice_cut_off", (DL_FUNC) &_rocmice_cut_off, 2},
    {"_rocmice_varlogit", (DL_FUNC) &_rocmice_varlogit, 2},
    {"_rocmice_apply_cut_off", (DL_FUNC) &_rocmice_apply_cut_off, 2},
    {"_rocmice_count_vals", (DL_FUNC) &_rocmice_count_vals, 2},
    {"_rocmice_apply_metrics", (DL_FUNC) &_rocmice_apply_metrics, 2},
    {"_rocmice_rcpp_rowMeans", (DL_FUNC) &_rocmice_rcpp_rowMeans, 1},
    {"_rocmice_rcpp_var", (DL_FUNC) &_rocmice_rcpp_var, 1},
    {"_rocmice_rowVars", (DL_FUNC) &_rocmice_rowVars, 1},
    {"_rocmice_pool_rr_m", (DL_FUNC) &_rocmice_pool_rr_m, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rocmice(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}