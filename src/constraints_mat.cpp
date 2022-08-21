#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <pairwise_diff_vec.cpp>
#include <cmath> 
using namespace Rcpp;
using namespace std;
using namespace Eigen; 
// [[Rcpp::depends(RcppEigen)]] 
// [[Rcpp::export]]
MatrixXd pairwise_mat_cpp(const MatrixXd &X, const int num_terms, const bool sym, VectorXi &sub_sample_indx){
  const int n(X.rows());
  const int p(X.cols());
  
  //Placeholders
  MatrixXd X_diff_mat(num_terms, p);
  
  //Compute pairwise difference each column of X
  for(int i = 0; i < p; ++i) {
    X_diff_mat.col(i) = pairwise_vec_cpp(X.col(i), num_terms, sym, sub_sample_indx);
  }
  return  X_diff_mat;
}

// [[Rcpp::export]]
Rcpp::List create_constr_mat(const MatrixXd &X, const int num_terms, const bool sym, VectorXi &sub_sample_indx){
  const int n(X.rows());
  const int p(X.cols());
  //Placeholders
  MatrixXd constr_mat(num_terms+2*p, 2*(num_terms+p));
  MatrixXd X_diff_mat = pairwise_mat_cpp(X, num_terms, sym, sub_sample_indx);
  
  //Create constraints matrix
  constr_mat <<  MatrixXd::Identity(num_terms, num_terms), -1*MatrixXd::Identity(num_terms, num_terms), MatrixXd::Zero(num_terms, p), X_diff_mat,
                 MatrixXd::Zero(p, 2*num_terms), MatrixXd::Identity(p, p), MatrixXd::Identity(p, p),
                 MatrixXd::Zero(p, 2*num_terms), MatrixXd::Identity(p, p), -1*MatrixXd::Identity(p, p);
  SparseMatrix<double> spMat = constr_mat.sparseView();
  return Rcpp::List::create(Rcpp::Named("constr_mat") = spMat,
                              Rcpp::Named("X_diff_mat") = X_diff_mat);
}

// [[Rcpp::export]]
Rcpp::List create_constr_mat_dense(const MatrixXd &X, const int num_terms, const bool sym, VectorXi &sub_sample_indx){
  const int n(X.rows());
  const int p(X.cols());
  //Placeholders
  MatrixXd constr_mat(num_terms+2*p, 2*(num_terms+p));
  MatrixXd X_diff_mat = pairwise_mat_cpp(X, num_terms, sym, sub_sample_indx);
  
  //Create constraints matrix
  constr_mat <<  MatrixXd::Identity(num_terms, num_terms), -1*MatrixXd::Identity(num_terms, num_terms), MatrixXd::Zero(num_terms, p), X_diff_mat,
                 MatrixXd::Zero(p, 2*num_terms), MatrixXd::Identity(p, p), MatrixXd::Identity(p, p),
                 MatrixXd::Zero(p, 2*num_terms), MatrixXd::Identity(p, p), -1*MatrixXd::Identity(p, p);
  return Rcpp::List::create(Rcpp::Named("constr_mat") = constr_mat,
                            Rcpp::Named("X_diff_mat") = X_diff_mat);
}

/*
// [[Rcpp::export]]
Rcpp::List create_constr_mat_ecos(const MatrixXd &X, const int num_terms, const bool sym, VectorXi &sub_sample_indx){
  const int n(X.rows());
  const int p(X.cols());
  //Placeholders
  MatrixXd constr_mat(num_terms, 2*(num_terms+p));
  MatrixXd X_diff_mat = pairwise_mat_cpp(X, num_terms, sym, sub_sample_indx);
  
  //Create constraints matrix
  constr_mat <<  MatrixXd::Identity(num_terms, num_terms), -1*MatrixXd::Identity(num_terms, num_terms), MatrixXd::Zero(num_terms, p), X_diff_mat;
  SparseMatrix<double> spMat = constr_mat.sparseView();
  return Rcpp::List::create(Rcpp::Named("constr_mat") = constr_mat,
                            Rcpp::Named("X_diff_mat") = X_diff_mat);
}

*/











