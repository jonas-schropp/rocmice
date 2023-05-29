// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// delongPlacementsCpp
List delongPlacementsCpp(List roc);
RcppExport SEXP _rocmice_delongPlacementsCpp(SEXP rocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type roc(rocSEXP);
    rcpp_result_gen = Rcpp::wrap(delongPlacementsCpp(roc));
    return rcpp_result_gen;
END_RCPP
}
// order_descending
std::vector<int> order_descending(std::vector<double> x);
RcppExport SEXP _rocmice_order_descending(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_descending(x));
    return rcpp_result_gen;
END_RCPP
}
// reorder
std::vector<int> reorder(const std::vector<int>& values, const std::vector<int>& indices);
RcppExport SEXP _rocmice_reorder(SEXP valuesSEXP, SEXP indicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type indices(indicesSEXP);
    rcpp_result_gen = Rcpp::wrap(reorder(values, indices));
    return rcpp_result_gen;
END_RCPP
}
// varproplogit
double varproplogit(double events, double n, double corr);
RcppExport SEXP _rocmice_varproplogit(SEXP eventsSEXP, SEXP nSEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type corr(corrSEXP);
    rcpp_result_gen = Rcpp::wrap(varproplogit(events, n, corr));
    return rcpp_result_gen;
END_RCPP
}
// logittrans
double logittrans(double x);
RcppExport SEXP _rocmice_logittrans(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logittrans(x));
    return rcpp_result_gen;
END_RCPP
}
// tpr
NumericVector tpr(IntegerVector r, IntegerVector p, double corr);
RcppExport SEXP _rocmice_tpr(SEXP rSEXP, SEXP pSEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type corr(corrSEXP);
    rcpp_result_gen = Rcpp::wrap(tpr(r, p, corr));
    return rcpp_result_gen;
END_RCPP
}
// fpr
NumericVector fpr(IntegerVector r, IntegerVector p, double corr);
RcppExport SEXP _rocmice_fpr(SEXP rSEXP, SEXP pSEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type corr(corrSEXP);
    rcpp_result_gen = Rcpp::wrap(fpr(r, p, corr));
    return rcpp_result_gen;
END_RCPP
}
// get_roc
NumericMatrix get_roc(std::vector<int> response, std::vector<double> prediction, double corr);
RcppExport SEXP _rocmice_get_roc(SEXP responseSEXP, SEXP predictionSEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type response(responseSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< double >::type corr(corrSEXP);
    rcpp_result_gen = Rcpp::wrap(get_roc(response, prediction, corr));
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
// combineVecs
std::vector<double> combineVecs(NumericVector A, NumericVector B);
RcppExport SEXP _rocmice_combineVecs(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(combineVecs(A, B));
    return rcpp_result_gen;
END_RCPP
}
// getUniqueValues
std::vector<double> getUniqueValues(std::vector<double> vec);
RcppExport SEXP _rocmice_getUniqueValues(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(getUniqueValues(vec));
    return rcpp_result_gen;
END_RCPP
}
// findIndex
int findIndex(double value, const NumericVector vec);
RcppExport SEXP _rocmice_findIndex(SEXP valueSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type value(valueSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(findIndex(value, vec));
    return rcpp_result_gen;
END_RCPP
}
// full_join_rcpp
NumericMatrix full_join_rcpp(NumericMatrix m, NumericVector fpri_vals, NumericVector zero);
RcppExport SEXP _rocmice_full_join_rcpp(SEXP mSEXP, SEXP fpri_valsSEXP, SEXP zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fpri_vals(fpri_valsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type zero(zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(full_join_rcpp(m, fpri_vals, zero));
    return rcpp_result_gen;
END_RCPP
}
// filterArray
IntegerVector filterArray(NumericMatrix mat);
RcppExport SEXP _rocmice_filterArray(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(filterArray(mat));
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
    {"_rocmice_delongPlacementsCpp", (DL_FUNC) &_rocmice_delongPlacementsCpp, 1},
    {"_rocmice_order_descending", (DL_FUNC) &_rocmice_order_descending, 1},
    {"_rocmice_reorder", (DL_FUNC) &_rocmice_reorder, 2},
    {"_rocmice_varproplogit", (DL_FUNC) &_rocmice_varproplogit, 3},
    {"_rocmice_logittrans", (DL_FUNC) &_rocmice_logittrans, 1},
    {"_rocmice_tpr", (DL_FUNC) &_rocmice_tpr, 3},
    {"_rocmice_fpr", (DL_FUNC) &_rocmice_fpr, 3},
    {"_rocmice_get_roc", (DL_FUNC) &_rocmice_get_roc, 3},
    {"_rocmice_varlogit", (DL_FUNC) &_rocmice_varlogit, 2},
    {"_rocmice_combineVecs", (DL_FUNC) &_rocmice_combineVecs, 2},
    {"_rocmice_getUniqueValues", (DL_FUNC) &_rocmice_getUniqueValues, 1},
    {"_rocmice_findIndex", (DL_FUNC) &_rocmice_findIndex, 2},
    {"_rocmice_full_join_rcpp", (DL_FUNC) &_rocmice_full_join_rcpp, 3},
    {"_rocmice_filterArray", (DL_FUNC) &_rocmice_filterArray, 1},
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
