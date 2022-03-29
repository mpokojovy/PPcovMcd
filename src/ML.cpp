//#include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]

inline static double log_double(double x){return ::log(x);}
inline static double sqr_double(double x) {return(x*x);}
inline static double sqrt_double(double x){return ::sqrt(x);}
inline static double pow_double(double base, double exp){return ::pow(base, exp);}

inline arma::vec pow2(const arma::vec &x) {
  return x % x;
}

// Zeroth derivative
inline arma::vec phi(const arma::vec &x, double mu, double sigma) {
  return arma::normpdf(x, mu, sigma);
};

// First derivatives
inline arma::vec phi_mu(const arma::vec &x, double mu, double sigma) {
  return ((x - mu)/sqr_double(sigma)) % arma::normpdf(x, mu, sigma);
};

inline arma::vec phi_sigma(const arma::vec &x, double mu, double sigma) {
  return ((pow2((x - mu)/sigma) - 1)/sigma) % arma::normpdf(x, mu, sigma);
};

// Second derivatives
inline arma::vec phi_mu_mu(const arma::vec &x, double mu, double sigma) {
  double sigma2 = sigma *sigma;
  return (pow2(x - mu) - sigma2)/(sigma2*sigma2) % arma::normpdf(x, mu, sigma);
};

inline arma::vec phi_sigma_sigma(const arma::vec &x, double mu, double sigma) {
  double sigma2 = sigma *sigma;
  double sigma4 = sigma2*sigma2;
  double sigma6 = sigma2*sigma4;

  arma::vec a = pow2(x - mu);
  return (pow2(a - sigma2)/sigma6 - (3*a - sigma2)/sigma4) % arma::normpdf(x, mu, sigma);
};

inline arma::vec phi_mu_sigma(const arma::vec &x, double mu, double sigma) {
  double sigma2 = sigma *sigma;
  double sigma3 = sigma2*sigma;

  arma::vec a = (x - mu);
  return (a/sigma3) % (pow2(a)/sigma2 - 3) % arma::normpdf(x, mu, sigma);
};

// Negative log-likelihood function
// [[Rcpp::export]]
Rcpp::List ML_objective(arma::vec par, const arma::vec &x, unsigned int nargout = 1) {
  unsigned int n = x.n_elem;

  double pi = par[0];
  double mu1 = par[1];
  double mu2 = par[2];
  double sigma1 = par[3];
  double sigma2 = par[4];

  double       obj = 0;
  arma::vec5  grad = arma::zeros<arma::vec>(5);
  arma::mat55 hess = arma::zeros<arma::mat>(5, 5);

  arma::vec phi1 = phi(x, mu1, sigma1);
  arma::vec phi2 = phi(x, mu2, sigma2);
  arma::vec f    = arma::zeros<arma::vec>(n);

  for (unsigned int i = 0; i < n; i++) {
    f[i] = (1 - pi)*phi1[i] + pi*phi2[i] + 1E-100;

    obj += log_double(f[i])/n;
  }

  if (nargout > 1) {
    arma::vec phi_mu1 = phi_mu(x, mu1, sigma1);
    arma::vec phi_mu2 = phi_mu(x, mu2, sigma2);
    arma::vec phi_sigma1 = phi_sigma(x, mu1, sigma1);
    arma::vec phi_sigma2 = phi_sigma(x, mu2, sigma2);

    grad[0] = arma::mean((-phi1 + phi2)/f);
    grad[1] = arma::mean((1 - pi)*phi_mu1/f);
    grad[2] = arma::mean((pi)*phi_mu2/f);
    grad[3] = arma::mean((1 - pi)*phi_sigma1/f);
    grad[4] = arma::mean((pi)*phi_sigma2/f);

    if (nargout > 2) {
      arma::vec phi_mu1_mu1       = phi_mu_mu(x, mu1, sigma1);
      arma::vec phi_mu1_sigma1    = phi_mu_sigma(x, mu1, sigma1);
      arma::vec phi_sigma1_sigma1 = phi_sigma_sigma(x, mu1, sigma1);

      arma::vec phi_mu2_mu2       = phi_mu_mu(x, mu2, sigma2);
      arma::vec phi_mu2_sigma2    = phi_mu_sigma(x, mu2, sigma2);
      arma::vec phi_sigma2_sigma2 = phi_sigma_sigma(x, mu2, sigma2);

      // partial pi pi
      hess(0, 0) = -arma::mean(pow2((-phi1 + phi2)/f));
      // partial pi partial mu1
      hess(0, 1) = -arma::mean((1 - pi)*(phi_mu1 % (-phi1 + phi2))/pow2(f) + phi_mu1/f);
      hess(1, 0) =  hess(0, 1);
      // partial pi partial mu2
      hess(0, 2) = -arma::mean(pi*(phi_mu2 % (-phi1 + phi2))/pow2(f) - phi_mu2/f);
      hess(2, 0) =  hess(0, 2);
      // partial pi partial sigma1
      hess(0, 3) = -arma::mean((1 - pi)*(phi_sigma1 % (-phi1 + phi2))/pow2(f) + phi_sigma1/f);
      hess(3, 0) =  hess(0, 3);
      // partial pi partial sigma2
      hess(0, 4) = -arma::mean(pi*(phi_sigma2 % (-phi1 + phi2))/pow2(f) - phi_sigma2/f);
      hess(4, 0) =  hess(0, 4);

      // partial mu1 partial mu1
      hess(1, 1) = -arma::mean(pow2((1 - pi)*phi_mu1/f) - (1 - pi)*phi_mu1_mu1/f);
      // partial mu2 partial mu2
      hess(2, 2) = -arma::mean(pow2(pi*phi_mu2/f) - pi*phi_mu2_mu2/f);
      // partial mu1 partial mu2
      hess(1, 2) = -arma::mean(pi*(1 - pi)*(phi_mu1 % phi_mu2)/pow2(f));
      // partial mu2 partial mu1
      hess(2, 1) = hess(1, 2);

      // partial sigma1 partial sigma1
      hess(3, 3) = -arma::mean(pow2((1 - pi)*phi_sigma1/f) - (1 - pi)*phi_sigma1_sigma1/f);
      // partial sigma2 partial sigma2
      hess(4, 4) = -arma::mean(pow2(pi*phi_sigma2/f) - pi*phi_sigma2_sigma2/f);
      // partial sigma1 partial sigma2
      hess(3, 4) = -arma::mean(pi*(1 - pi)*(phi_sigma1 % phi_sigma2)/pow2(f));
      // partial sigma2 partial sigma1
      hess(4, 3) = hess(3, 4);

      // partial mu1 partial sigma1
      hess(1, 3) = -arma::mean(pow2((1 - pi)/f) % phi_mu1 % phi_sigma1 - (1 - pi)*phi_mu1_sigma1/f);
      hess(3, 1) =  hess(1, 3);
      // partial mu2 partial sigma2
      hess(2, 4) = -arma::mean(pow2(pi/f) % phi_mu2 % phi_sigma2 - pi*phi_mu2_sigma2/f);
      hess(4, 2) =  hess(2, 4);
      // partial mu1 partial sigma2
      hess(1, 4) = -arma::mean(pi*(1 - pi)*(phi_mu1 % phi_sigma2)/pow2(f));
      hess(4, 1) =  hess(1, 4);
      // partial mu2 partial sigma1
      hess(2, 3) = -arma::mean(pi*(1 - pi)*(phi_mu2 % phi_sigma1)/pow2(f));
      hess(3, 2) =  hess(2, 3);
    }
  }

  if (nargout > 3) {
    arma::vec newton_dir; // Newton descent direction

    double hess_rcond = arma::rcond(hess);

    if (std::isnan(hess_rcond)) {
      newton_dir = -grad;
    } else if (hess_rcond >= 1E-9) {
      newton_dir = -arma::solve(hess, grad);
    } else {
      newton_dir = -grad;
    }

    return Rcpp::List::create(Rcpp::Named("obj")  = -obj,
                              Rcpp::Named("grad") = -grad,
                              Rcpp::Named("hess") = -hess,
                              Rcpp::Named("newton.dir") = -newton_dir);
  } else if (nargout > 2) {
    return Rcpp::List::create(Rcpp::Named("obj")  = -obj,
                              Rcpp::Named("grad") = -grad,
                              Rcpp::Named("hess") = -hess);
  } else if (nargout > 1) {
    return Rcpp::List::create(Rcpp::Named("obj")  = -obj,
                              Rcpp::Named("grad") = -grad);
  } else {
    return Rcpp::List::create(Rcpp::Named("obj")  = -obj);
  }
}

// Projector
// [[Rcpp::export]]
arma::vec5 ML_projector(arma::vec par, double sigma_min)
{
  par[0] = std::max(0.0, std::min(par[0], 1.0));
  par[3] = std::max(par[3], sigma_min);
  par[4] = std::max(par[4], sigma_min);

  return par;
}
