/* -*- mode: C++
 * 
 * Computes the 1D/univariate Rousseeuw's MCD estimator and its gradient as a function of the sample
 * Reference: Rousseeuw, P. & Leroy, A. (1987). Robust regression and outlier detection. Wiley, New York, NY, pp. 171-172
 * The "gradient" is computed under assumption the samples contains *no* ties.
 * In case of ties, the output can be viewed as a non-trivial element of the semigradient, not the gradient.
 * 
 * (C) Michael Pokojovy, 2019  
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//#define M_PI        3.141592653589793238462643383280 

inline static double sqr_double(double x) {return(x*x);}
inline static double log_double(double x){return ::log(x);}
inline static double sqrt_double(double x){return ::sqrt(x);}
inline static double abs_double(double x){return std::abs(x);}

// [[Rcpp::export]]
Rcpp::List mcd1D(const arma::vec &x, unsigned int h) {
  const int n = x.n_elem;

  if (h == n) {
    arma::uvec rank = arma::stable_sort_index(x); // order low to high

    return Rcpp::List::create(Rcpp::Named("mean") = mean(x), Rcpp::Named("var") = arma::var(x),
                              Rcpp::Named("best") = rank,
                              Rcpp::Named("h") = h);
  } else {
    arma::uvec rank = arma::stable_sort_index(x); // order low to high
    arma::vec y = x(rank); // sorted sample
    
    arma::vec y_bar(n - h + 1);
    arma::vec s_sq(n - h + 1);

    y_bar(0) = arma::mean(y.subvec(0, h - 1));
    s_sq(0)  = arma::var(y.subvec(0, h - 1))*(h - 1);

    for (unsigned int j = 1; j < n - h + 1; j++) {
      y_bar(j) = (h*y_bar(j - 1) - y(j - 1) + y(j + h - 1))/h;
      s_sq(j)   = s_sq(j - 1) - sqr_double(y(j - 1)) + sqr_double(y(j + h - 1)) - h*sqr_double(y_bar(j)) + h*sqr_double(y_bar(j - 1));
    }

    const unsigned int j_ast = s_sq.index_min();

    return Rcpp::List::create(Rcpp::Named("mean") = y_bar(j_ast), Rcpp::Named("var") = s_sq(j_ast)/(h - 1),
                              Rcpp::Named("best") = rank.subvec(j_ast, j_ast + (h - 1)) + 1,
                              Rcpp::Named("h") = h);
  }
}

// [[Rcpp::export]]
arma::vec mcd1D_grad(const arma::vec &x, const arma::uvec &best) {
  const int n = x.n_elem;
  const int h = best.n_elem;
  
  const double mean = arma::mean(x(best - 1));

  arma::vec grad(n);
  grad.zeros();
  
  grad(best - 1) = (x(best - 1) - mean)/(0.5*(h - 1));

  return grad;
}

// [[Rcpp::export]]
double PP_objective(const arma::vec &xp, const arma::uvec &best) {
  const int n = xp.n_elem;
  
  const double s2 = arma::var(xp);

  // TODO: Implement a better smallness condition
  if (s2 == 0.0)
    return R_PosInf;
  else {
    const double mcd_s2 = arma::var(xp(best - 1));
    
    // TODO: Implement a better smallness condition
    if (mcd_s2 == 0.0)
      return R_PosInf;
    else
      return log_double(s2/mcd_s2);
  }
}

// [[Rcpp::export]]
arma::colvec PP_gradient(const arma::vec &xp, const arma::mat &x, const arma::uvec &best) {
  const int n = xp.n_elem;
  const int k = x.n_cols;
  const int h = best.n_elem;

  double s2 = arma::var(xp);

  // TODO: Implement a better smallness condition
  if (s2 == 0.0) {
    return arma::zeros<arma::rowvec>(k);
  } else {
    const double xbar = arma::mean(xp);
    
    const double mcd_xbar = arma::mean(xp(best - 1));
    const double mcd_s2   = arma::var(xp(best - 1));

    // TODO: Implement a better smallness condition
    if (mcd_s2 == 0)
      return arma::zeros<arma::rowvec>(k);
    else {
      arma::vec grad_var = (xp - xbar)/(0.5*(n - 1));
      arma::vec grad_mcd = mcd1D_grad(xp, best);
      
      arma::mat grad = (arma::mat) (grad_var/s2 - grad_mcd/mcd_s2);

      return (arma::colvec) (grad.t()*x).t();
    }
  }
}

// [[Rcpp::export]]
arma::colvec PP_project(const arma::colvec &dir) {
  return dir/sqrt_double(arma::dot(dir, dir)); 
}
