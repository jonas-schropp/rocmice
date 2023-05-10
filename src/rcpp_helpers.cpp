#include <Rcpp.h>
using namespace Rcpp;

//' TPR
//' @param r response vector, 0 / 1
//' @param p prediction vector, 0 / 1
//' 
//' @keywords Internal
// [[Rcpp::export]]
double tpr(IntegerVector r, IntegerVector p) {
  int n = r.size();
  int num_positives = 0;
  int true_positives = 0;

  for (int i = 0; i < n; i++) {
    if (r[i] == 1) {
      if (p[i] == 1) {
        true_positives++;
      }
      num_positives++;
    }
  }

  return (double) true_positives / num_positives;
}

//' FPR
//' @param r response vector, 0 / 1
//' @param p prediction vector, 0 / 1
//' 
//' @keywords Internal
// [[Rcpp::export]]
double fpr(IntegerVector r, IntegerVector p) {
  int n = r.size();
  int num_negatives = 0;
  int false_positives = 0;

  for (int i = 0; i < n; i++) {
    if (r[i] == 0) {
      if (p[i] == 1) {
        false_positives++;
      }
      num_negatives++;
    }
  }

  return (double) false_positives / num_negatives;
}

//' Variance for a proportion
//' @param est vector of tpr or fpr
//' @param n number of positives in r or negatives in r
//' 
//' @keywords Internal
// [[Rcpp::export]]
NumericVector var_prop(NumericVector est, double n) {
  int len = est.length();
  NumericVector variances(len);

  for (int i = 0; i < len; i++) {
    variances[i] = (est[i] * (1 - est[i])) / n;
  }

  return variances;
}

//' check if a value is bigger than a cutoff
//' @param score NumericVector of scores
//' @param cutoff double, the cutoff
//' 
//' @keywords Internal
// [[Rcpp::export]]
IntegerVector cut_off(NumericVector score, double cutoff) {
  int len = score.length();
  IntegerVector result(len);

  for (int i = 0; i < len; i++) {
    result[i] = (score[i] > cutoff) ? 1 : 0;
  }

  return result;
}

//' transform variances to logit scale
//' @param var NumericVector of variances
//' @param est NumericVector of proportions
//' 
//' @keywords Internal
// [[Rcpp::export]]
NumericVector varlogit(NumericVector var, NumericVector est) {
  int len = var.length();
  NumericVector res(len);

  for (int i = 0; i < len; i++) {
    if (est[i] == 0) {
      res[i] = 1e-8;
    } else {
      res[i] = sqrt(var[i] / (est[i] * (1 - est[i])));
    }
  }

  return res;
}


//' Apply all possible cutoffs 
//' 
//' @param score the score
//' @param unique_vals the possible cutoffs
//' 
//' @return 
//' A matrix of dimensions length(score) x length(cutoffs)
//' 
//' 
//' @keywords Internal
// [[Rcpp::export]]
NumericMatrix apply_cut_off(NumericVector score, NumericVector unique_vals) {
  int len = unique_vals.length() - 1;
  int score_len = score.length();
  NumericMatrix res(score_len, len);

  for (int i = 0; i < len; i++) {
    res(_, i) = cut_off(score, unique_vals[i]);
  }

  return res;
}

//' Function to count number of val in x
//' @param x IntegerVector
//' @param val int, the value to count
//' 
//' @keywords Internal
// [[Rcpp::export]]
int count_vals(IntegerVector x, int val) {
  int count = 0;
  for (int i = 0; i < x.length(); i++) {
    if (x[i] == val) {
      count++;
    }
  }
  return count;
}


//' Function to apply metrics
//'
//' @param m matrix with one column for each cutoff and one row for each 
//' value in score
//' @param r binary response vector
//' 
//' @return 
//' a list with 3 elements:
//' \item{tpri}{TPR values}
//' \item{fpri}{FPR values}
//' \item{vartpri}{variance for TPR}
//' 
//' @keywords Internal
//'
// [[Rcpp::export]]
List apply_metrics(const IntegerMatrix& m, const IntegerVector& r) {
  int len = m.ncol();
  NumericVector tpri(len);
  NumericVector fpri(len);
  NumericVector vartpri(len);
  
  int pos = count_vals(r, 1);
  
  for (int i = 0; i < len; i++) {
    IntegerVector col = m(_, i);
    tpri[i] = tpr(r, col);
    fpri[i] = fpr(r, col);
  }
  
  vartpri = var_prop(tpri, pos);
  
  return List::create(
    Named("tpri") = tpri,
    Named("fpri") = fpri,
    Named("vartpri") = vartpri
  );
}


//' Equivalent to R rowMeans
//' @param x a matrix
//' @keywords Internal
// [[Rcpp::export]]
NumericVector rcpp_rowMeans(NumericMatrix x) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  NumericVector out(nrow);
  for(int i=0; i<nrow; i++) {
    double sum = 0;
    for(int j=0; j<ncol; j++) {
      sum += x(i, j);
    }
    out[i] = sum / ncol;
  }
  return out;
}


//' Equivalent to R var
//' @param samples a vector
//' @keywords Internal
// [[Rcpp::export]]
double rcpp_var(NumericVector samples)
{
     int size = samples.size();

     double variance = 0;
     double t = samples[0];
     for (int i = 1; i < size; i++)
     {
          t += samples[i];
          double diff = ((i + 1) * samples[i]) - t;
          variance += (diff * diff) / ((i + 1.0) *i);
     }

     return variance / (size - 1);
}


//' Apply var to matrix rows
//' @param x a matrix
//' @keywords Internal
// [[Rcpp::export]]
NumericVector rowVars(NumericMatrix x) {

  int nrow = x.nrow();
  NumericVector out(nrow);

  for(int i=0; i<nrow; i++) {
    NumericVector vec = x(i,_);
    out[i] = rcpp_var(vec);
  }

  return(out);

}




//' Pool mean and variance for many rows following Rubin's Rules
//' @param est a matrix of (logit-transformed) proportion by imp
//' @param var a matrix of (logit transformed) variance by imp
//' 
//' @keywords Internal
// [[Rcpp::export]]
List pool_rr_m (NumericMatrix est, NumericMatrix var) {

  int m = est.ncol();
  int l = est.nrow();
  NumericVector var_T(l);
  NumericVector var_w(l);
  NumericVector var_b(l);
  NumericVector mean_est(l);

  mean_est = rcpp_rowMeans(est);
  var_w = rcpp_rowMeans(var);
  var_b = rowVars(est);
  //var_b[is.nan(var_b)] <- 0

  for (int i=0; i<l; i++) {
    var_T[i] = var_w[i] + (1 + 1.0 / m) * var_b[i];
  }
  
  return List::create(
    Named("roc") = mean_est,
    Named("var_roc") = var_T
  );

}