% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pool_roc_rr.R
\name{interpol}
\alias{interpol}
\title{Interpolate values for tpri}
\usage{
interpol(l, fpr_vals, zero)
}
\arguments{
\item{l}{a matrix with 4 columns: "fpri", "tpri", "varfpri", "vartpri" as 
returned by `apply_metrics`}

\item{fpr_vals}{the FPR values to interpolate for}

\item{zero}{values to use for 0 cells}
}
\value{
a data.frame with 3 elements:
\item{tpri}{logit-transformed TPR values}
\item{fpri}{logit-transformed FPR values passed on from `fpr_vals`}
\item{varfpri}{logit-transformed variance for FPR}
\item{vartpri}{logit-transformed variance for TPR}
}
\description{
Unlike the method for meta-analysis by Martínez-Camblor, P. (2016), this 
function does not use linear interpolation because we have enough values 
for tpr available to simply use the values in the data.
}
\keyword{Internal}
