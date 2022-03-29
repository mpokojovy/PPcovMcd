#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

inline static double pow_double(double base, double exp){return ::pow(base, exp);}
inline static double sqrt_double(double x) {return ::sqrt(x);}

// [[Rcpp::export]]
arma::mat arma_ast(const arma::mat &A, const arma::mat &B) {
  return A*B;
}

// [[Rcpp::export]]
Rcpp::List mah_standard(const arma::mat &x, const arma::colvec &loc, const arma::mat &cov, unsigned int rank) {
  const unsigned int n = x.n_rows, p = x.n_cols;

  arma::vec eigval;
  arma::mat U, V;

  arma::svd(U, eigval, V, cov);

  bool singular = false;

  // root.cov
  arma::mat root_cov = U;

  arma::vec d = arma::zeros<arma::colvec>(p);
  for (unsigned int i = 0; i < rank; i++) {
    if (eigval[p - 1 - i] > 0)
      d[p - 1 - i] = sqrt_double(eigval[p - 1 - i]);
    else
      singular = true;
  }

  for (unsigned int i = 0; i < p; i++) {
    root_cov.col(i) *= d[i];
  }

  root_cov = root_cov*V.t();

  // inv.root.cov
  arma::mat inv_root_cov = U;

  for (unsigned int i = 0; i < rank; i++) {
    if (eigval[p - 1 - i] > 0)
      d[p - 1 - i] = 1/d[p - 1 - i];
    else
      singular = true;
  }

  for (unsigned int i = 0; i < p; i++) {
    inv_root_cov.col(i) *= d[i];
  }

  inv_root_cov = inv_root_cov*V.t();

  // z.scores
  arma::mat z_scores = arma::trans(x);
  for (unsigned int i = 0; i < n; i++) {
    z_scores.col(i) -= loc;
  }

  z_scores = inv_root_cov*z_scores;

  // return
  return Rcpp::List::create(Rcpp::Named("z.scores")     = arma::trans(z_scores),
                            Rcpp::Named("loc")          = loc,
                            Rcpp::Named("root.cov")     = root_cov,
                            Rcpp::Named("inv.root.cov") = inv_root_cov,
                            Rcpp::Named("singular")     = singular);
}

// [[Rcpp::export]]
arma::mat proj_orth_comp(const arma::mat &x, const arma::colvec &dir) {
  const unsigned int n = x.n_rows;
  
  arma::mat y = arma::trans(x);
  
  for (unsigned int i = 0; i < n; i++) {
    y.col(i) -= arma::cdot(y.col(i), dir)*dir;
  }
  
  return arma::trans(y);
}
