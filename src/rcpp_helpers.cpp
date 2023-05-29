#include <Rcpp.h>
using namespace Rcpp;


//' C++ implementation of R `order`
//' Assuming we don't have any missing values, only one vector to order and 
//' always want to order descending.
//' @param x double vector to order
//' @keywords Internal
// [[Rcpp::export]]
std::vector<int> order_descending(std::vector<double> x) {
  int n = x.size();

  // Create an index vector
  std::vector<int> indices(n);
  for (int i = 0; i < n; ++i) {
    indices[i] = i;
  }

  // Sort the indices based on the values in x in descending order
  std::sort(indices.begin(), indices.end(), [&x](int i, int j) {
    return x[i] > x[j];
  });

  return indices;
}


//' Reorder a int vector by indices
//' equivalent to R vector[indices]
//' @param values the vector to reorder/subset
//' @param indices integer vector of indices to order by
// [[Rcpp::export]]
std::vector<int> reorder(const std::vector<int>& values, const std::vector<int>& indices) {

  int n = indices.size();
  std::vector<int> result(n);

  for (int i = 0; i < n; ++i) {
    result[i] = values[indices[i]];  
  }

  return result;
}


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

//' Logit transformation of a double
//' @param x double to transform
// [[Rcpp::export]]
double logittrans(double x) {

  double res = std::log(x / (1 - x));
  return res;

}


//' TPR
//' @param r response vector, 0 / 1
//' @param p prediction vector, 0 / 1
//' @param corr continuity correction
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



//' Calculates ROC curve (tpr and fpr at each cutoff)
//'
//' Using a fast greedy algorithm
//'
//' @param response integer, response vector
//' @param prediction double, predictions
//' @param corr Continuity correction
//' 
//' @keywords Internal
// [[Rcpp::export]]
NumericMatrix get_roc (
    std::vector<int> response, 
    std::vector<double> prediction, 
    double corr
    ) {

  // initialize stuff
  int n = response.size();
  //double tmp = 1.0;
  NumericMatrix result(n + 1, 4);

  // tpr calculation is faster if we order response first
  std::vector<int> index = order_descending(prediction);
  std::vector<int> ordered_response = reorder(response, index);

  //std::vector<int> tpc = cumsum_rcpp(ordered_response, corr);
  double events = std::accumulate(response.begin(), response.end(), 0) * 1.0;
  double non_events = n - events;

  events = events + 2 * corr;
  non_events = non_events + 2 * corr;

  // initialize vectors to fill
  NumericVector tpc(n + 1);
  NumericVector tpr(n + 1);
  NumericVector fpc(n + 1);
  NumericVector fpr(n + 1);
  NumericVector var_tpr(n + 1);
  NumericVector var_fpr(n + 1);

  // add 0 at the start so curve always starts from 0
  tpc[0] = 0.0;
  fpc[0] = 0.0;
  tpr[0] = logittrans(corr / events);
  fpr[0] = logittrans(corr / non_events);
  var_tpr[0] = varproplogit(0.0, events, corr);
  var_fpr[0] = varproplogit(0.0, non_events, corr);

  // Main part - loop to fill values:
  for (int i = 1; i < tpc.size(); i++) {

    tpc[i] = tpc[i - 1] + ordered_response[i - 1];
    fpc[i] = fpc[i - 1] - (ordered_response[i - 1] - 1);

    // logit transform
    tpr[i] = logittrans((tpc[i] + corr) / events);
    fpr[i] = logittrans((fpc[i] + corr) / non_events);

    // Already transformed variance of TPR and FPR
    var_tpr[i] = varproplogit(tpc[i], events, corr);
    var_fpr[i] = varproplogit(fpc[i], non_events, corr);

  }

  // Return as matrix with 4 columns
  result(_, 0) = fpr;
  result(_, 1) = tpr;
  result(_, 2) = var_fpr;
  result(_, 3) = var_tpr;

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



