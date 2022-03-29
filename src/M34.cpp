//#define ARMA_DONT_USE_WRAPPER
//#include <Rcpp.h>
#include <RcppArmadillo.h>
//#include <cmath>
//#include <math.h>
//#include <Rmath.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

inline static double sqrt_double(double x) {return ::sqrt(x);}

// [[Rcpp::export]]
Rcpp::List M3_vec(const arma::mat &x, const arma::rowvec &xbar, const arma::mat &S) {
  const unsigned int n = x.n_rows, p = x.n_cols;
  
  arma::rowvec v = arma::zeros<arma::rowvec>(p);
  
  arma::mat S_inv = arma::inv(S);
  
  for (unsigned int i = 0; i < n; i++) {
    arma::rowvec y = x.row(i) - xbar;
    double d2 = arma::dot(y*S_inv, y);
    v += (d2/n)*y;
  }
  
  v = v/sqrt_double(arma::dot(v, v));
  return Rcpp::List::create(Rcpp::Named("v1") = v);
}

// [[Rcpp::export]]
arma::mat M4_matrix(const arma::mat &x, const arma::rowvec &xbar, const arma::mat &S) {
  const unsigned int n = x.n_rows, p = x.n_cols;
  
  arma::mat M = arma::zeros<arma::mat>(p, p);
  
  arma::mat S_inv = arma::inv(S);
  
  for (unsigned int i = 0; i < n; i++) {
    arma::rowvec y = x.row(i) - xbar;
    double d2 = arma::dot(y*S_inv, y);
    M += (d2/n)*y.t()*y;
  }
  
  M = 0.5*(M + M.t());
  return M;
}

// [[Rcpp::export]]
Rcpp::List M4_vec(const arma::mat &x, const arma::rowvec &xbar, const arma::mat &S) {
  arma::mat M = M4_matrix(x, xbar, S);

  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, M);

  unsigned int i1 = eigval.index_min();
  unsigned int ip = eigval.index_max();

  return Rcpp::List::create(Rcpp::Named("v1") = eigvec.col(i1), Rcpp::Named("vp") = eigvec.col(ip));
}

// [[Rcpp::export]]
Rcpp::List M34_vec(const arma::mat &x, const arma::rowvec &xbar, const arma::mat &S) {
  const unsigned int n = x.n_rows, p = x.n_cols;
  
  arma::mat S_inv = arma::inv(S);
  
  arma::rowvec v = arma::zeros<arma::rowvec>(p);
  arma::mat M = arma::zeros<arma::mat>(p, p);
  
  for (unsigned int i = 0; i < n; i++) {
    arma::rowvec y = x.row(i) - xbar;
    double d2 = arma::dot(y*S_inv, y);
    v += (d2/n)*y;
    M += (d2/n)*y.t()*y;
  }
  
  v = v/sqrt_double(arma::dot(v, v));
  M = 0.5*(M + M.t());
  
  arma::vec eigval;
  arma::mat eigvec;
  
  arma::eig_sym(eigval, eigvec, M);
  
  unsigned int i1 = eigval.index_min();
  unsigned int ip = eigval.index_max();
  
  return Rcpp::List::create(Rcpp::Named("v1") = v,
                            Rcpp::Named("v2") = eigvec.col(i1), 
                            Rcpp::Named("v3") = eigvec.col(ip));
}
