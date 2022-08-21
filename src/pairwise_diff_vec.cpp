#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
using namespace Rcpp;
using namespace std;
using namespace Eigen; 
// [[Rcpp::depends(RcppEigen)]] 

// [[Rcpp::export]]
VectorXd pairwise_vec_cpp(const VectorXd &v, const int num_terms, const bool sym, const VectorXi &sub_sample_indx){
  const int n(v.size());
  int num_terms_all;
  if (sym==false && num_terms < n*(n-1)/2) {
    num_terms_all = n*(n-1)/2;
  } else {
    num_terms_all = num_terms;
  }
  
  //Placeholders
  VectorXd X_diff_vec(num_terms_all);
  int v_indx = 0;

  //Compute pariwise differences
  for(int i = 0; i < n; ++i) {
    if (sym == TRUE) {
      for(int j = 0; j < n; ++j) {
        if (i==j) {
          continue;
        } else {
          X_diff_vec(v_indx) = v(i) - v(j);
          v_indx++;
        }
      }
    } else {
      for(int j = i; j < n; ++j) {
        if (i==j) {
          continue;
        } else {
          X_diff_vec(v_indx) = v(i) - v(j);
          v_indx++;
        }
      }
    }
  }
  if (num_terms_all == num_terms) {
    return X_diff_vec;
  } else {
    VectorXd X_diff_vec_sub(num_terms);
    for (int i = 0; i < num_terms; ++i) {
      X_diff_vec_sub(i) = X_diff_vec(sub_sample_indx(i));
    }
    return X_diff_vec_sub;
  }
}

