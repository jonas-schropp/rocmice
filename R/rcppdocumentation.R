#' TPR
#' @name tpr
#' @param r response vector, 0 / 1
#' @param p prediction vector, 0 / 1
#' 
#' @keywords Internal
tpr

#' FPR
#' @name fpr
#' @param r response vector, 0 / 1
#' @param p prediction vector, 0 / 1
#' 
#' @keywords Internal
fpr

#' Variance for a proportion
#' @name var_prop
#' @param est vector of tpr or fpr
#' @param n number of positives in r or negatives in r
#' 
#' @keywords Internal
var_prop

#' check if a value is bigger than a cutoff
#' @name cut_off
#' @param score NumericVector of scores
#' @param cutoff double, the cutoff
#' 
#' @keywords Internal
cut_off

#' transform variances to logit scale
#' @name varlogit
#' @param var NumericVector of variances
#' @param est NumericVector of proportions
#' 
#' @keywords Internal
varlogit


#' Apply all possible cutoffs 
#' @name apply_cut_off
#' 
#' @param score the score
#' @param unique_vals the possible cutoffs
#' 
#' @return 
#' A matrix of dimensions length(score) x length(cutoffs)
#' 
#' 
#' @keywords Internal
apply_cut_off

#' Function to count number of val in x
#' @name count_vals
#' @param x IntegerVector
#' @param val int, the value to count
#' 
#' @keywords Internal
count_vals


#' Function to apply metrics
#' @name apply_metrics
#'
#' @param m matrix with one column for each cutoff and one row for each 
#' value in score
#' @param r binary response vector
#' 
#' @return 
#' a list with 3 elements:
  #' \item{tpri}{TPR values}
#' \item{fpri}{FPR values}
#' \item{vartpri}{variance for TPR}
#' 
#' @keywords Internal
#'
apply_metrics

#' Equivalent to R rowMeans
#' @name rcpp_rowMeans
#' @param x a matrix
#' @keywords Internal
rcpp_rowMeans

#' Equivalent to R var
#' @name rcpp_var
#' @param samples a vector
#' @keywords Internal
rcpp_var


#' Apply var to matrix rows
#' @name rowVars
#' @param x a matrix
#' @keywords Internal
rowVars




#' Pool mean and variance for many rows following Rubin's Rules
#' @name pool_rr_m
#' @param est a matrix of (logit-transformed) proportion by imp
#' @param var a matrix of (logit transformed) variance by imp
#' 
#' @keywords Internal
pool_rr_m