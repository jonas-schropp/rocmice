#include <Rcpp.h>
using namespace Rcpp;


//' Variance for a logit-transformed proportion
//' @param events double, number of events
//' @param n double, number of trials
//' @param corr continuity correction
//' 
//' @keywords Internal
// [[Rcpp::export]]
double varproplogit(double events, double n, double corr) {

  double res = 1 / (events + corr) + 1 / ((n - events + corr));
  return res;

}


// [[Rcpp::export]]
double logittrans(double x) {

  double res = std::log(x / (1 - x));
  return res;

}


//' TPR
//' @param r response vector, 0 / 1
//' @param p prediction vector, 0 / 1
//' 
//' @keywords Internal
// [[Rcpp::export]]
NumericVector tpr(IntegerVector r, IntegerVector p, double corr) {

  NumericVector res(2);
  int n = r.size();
  int num_positives = 0.0;
  int true_positives = 0.0;

  for (int i = 0; i < n; i++) {
    if (r[i] == 1) {
      if (p[i] == 1) {
        true_positives++;
      }
      num_positives++;
    }
  }
  
  double TPR = (true_positives + corr) / (num_positives + 2*corr);
  res[0] = logittrans(TPR);
  res[1] = varproplogit(true_positives, num_positives, corr);

  return res;
}

//' FPR
//' @param r response vector, 0 / 1
//' @param p prediction vector, 0 / 1
//' @param corr continuity correction
//' 
//' @keywords Internal
// [[Rcpp::export]]
NumericVector fpr(IntegerVector r, IntegerVector p, double corr) {
  
  NumericVector res(2);
  int n = r.size();
  int num_negatives = 0.0;
  int false_positives = 0.0;

  for (int i = 0; i < n; i++) {
    if (r[i] == 0) {
      if (p[i] == 1) {
        false_positives++;
      }
      num_negatives++;
    }
  }
  
  double FPR = (false_positives + corr) / (num_negatives + 2*corr);
  res[0] = logittrans(FPR);
  res[1] = varproplogit(false_positives, num_negatives, corr);

  return res;
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
    res[i] = (1 / (est[i] * (1 - est[i]))) * (1 / (est[i] * (1 - est[i]))) * var[i];
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
NumericMatrix apply_metrics(const IntegerMatrix& m, const IntegerVector& r, double corr) {

  int len = m.ncol();
  NumericVector tpri(len);
  NumericVector fpri(len);
  NumericVector vartpri(len);
  NumericVector varfpri(len);
  
  //int pos = count_vals(r, 1);
  
  for (int i = 0; i < len; i++) {
    IntegerVector col = m(_, i);
    NumericVector tmp0 = tpr(r, col, corr);
    NumericVector tmp1 = fpr(r, col, corr);
    tpri[i] = tmp0[0];
    fpri[i] = tmp1[0];
    vartpri[i] = tmp0[1];
    varfpri[i] = tmp1[1];
  }
  
  NumericMatrix result(len, 4);
  result(_, 0) = fpri;
  result(_, 1) = tpri;
  result(_, 2) = varfpri;
  result(_, 3) = vartpri;

  return result;
}


// Function to combine two NumericVectors
// [[Rcpp::export]]
std::vector<double> combineVecs(NumericVector A, NumericVector B) {
  
  std::vector<double> res(A.size() + B.size()); // preallocate memory
  
  res.insert( res.end(), A.begin(), A.end() );
  res.insert( res.end(), B.begin(), B.end() );
  
  
  return res;
  
}


// Function to get unique values from two NumericVectors
// [[Rcpp::export]]
std::vector<double> getUniqueValues(std::vector<double> vec) {
  
  std::vector<double>::iterator ip;
  int l = vec.size();
  
  // Using std::unique
  ip = std::unique(vec.begin(), vec.begin() + l);
  // * means undefined
  
  // Resizing the vector so as to remove the undefined terms
  vec.resize(std::distance(vec.begin(), ip));
  
  return vec;
}


// Function to sort a NumericVector in ascending order
//std::vector<double> sortVectorAscending(std::vector<double> vec) {
  
//  std::sort(vec.begin(), vec.end());
  
//  return vec;
//}


// Function to find the index of a double value in a NumericVector
// Returns NA if the value is not found
// [[Rcpp::export]]
int findIndex(double value, const NumericVector vec) {
  for (int i = 0; i < vec.size(); i++) {
    if (vec[i] == value) {
      return i;
    } 
  }
  int x = -999;
  return x;
}


// C++ function to create the matrix
// [[Rcpp::export]]
NumericMatrix full_join_rcpp(
    NumericMatrix m, 
    NumericVector fpri_vals, 
    NumericVector zero
) {
  
  // Get unique values of fpri_vals and fpri in m
  std::vector<double> fpri_all = combineVecs(m(_,0), fpri_vals);
  std::vector<double> fpri = getUniqueValues(fpri_all);
  
  NumericMatrix result_matrix(fpri.size(), 4);
  
  //std::vector<double> fpri_sorted = sortVectorAscending(fpri);
  // Sort the vector in ascending order
  std::sort(fpri.begin(), fpri.end());
  
  NumericVector fpriNV(fpri.begin(), fpri.end()); // Convert to NumericVector
  result_matrix(_, 0) = fpriNV;

  
  // Filling columns 2, 3, and 4 with corresponding values
  for (int i = 0; i < fpriNV.size(); i++) {
    double fpr = result_matrix(i, 0);
    int index = findIndex(fpr, m(_,0));
    
    if (index >= 0) {
      result_matrix(i, 1) = m(index, 1);
      result_matrix(i, 2) = m(index, 2);
      result_matrix(i, 3) = m(index, 3);
    } else if (index < 0 && i > 0) {
      result_matrix(i, 1) = result_matrix(i-1, 1);
      result_matrix(i, 2) = result_matrix(i-1, 2);
      result_matrix(i, 3) = result_matrix(i-1, 3);
    } else if (index < 0 && i == 0) {
      result_matrix(i, 1) = zero[0];
      result_matrix(i, 2) = zero[1];
      result_matrix(i, 3) = zero[2];
    }
    
  }
  
  return result_matrix;
}


// Function to remove duplicate values and retain only the max tpr for each fpr
// [[Rcpp::export]]
IntegerVector filterArray(NumericMatrix mat) {
  
  // Unfortunately have to clone to avoid side effects...
  NumericMatrix m = clone(mat);
  IntegerVector res(m.nrow());
  
  for (int i = 1; i < m.nrow(); i++) {

    
    if (m(i,0) != m(i-1,0)) {
      continue;
    } else if (m(i,0) == m(i-1,0)){
      if (m(i,1) >= m(i-1,1)) {
        res[i-1] = 1;
      } else if (m(i,1) < m(i-1,1)) {
        NumericVector tmp = m(i,_);
        m(i,_) = m(i-1,_);
        m(i-1,_) = tmp;
        res[i] = 1;
      } 
    } 
  }
  
  return res;
  
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
//' @param samples NumericVector
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



