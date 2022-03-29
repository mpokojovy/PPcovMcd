ML_bulk <- function(x, h)
{
  init.theta <- function(x1, x2) {
    n1    = length(x1); n2 = length(x2)

    mu    = mean(c(x1, x2))
    sigma = sd(c(x1, x2))

    mu1 = mean(x1)
    mu1 = if (is.na(mu1)) mu else mu1

    mu2 = mean(x2)
    mu2 = if (is.na(mu2)) mu else mu2

    sigma1 = sd(x1)
    sigma1 = if (is.na(sigma1)) sigma else if (sigma1 == 0.0) sigma else sigma1

    sigma2 = sd(x2)
    sigma2 = if (is.na(sigma2)) sigma else if (sigma2 == 0.0) sigma else sigma2

    if (n1 >= n2)
      return(c(n2/(n1 + n2), c(mu1, mu2), c(sigma1, sigma2)))
    else
      return(c(n1/(n1 + n2), c(mu2, mu1), c(sigma2, sigma1)))
  }

  xbar   = mean(x)
  sqrts2 = sd(x)

  mcd = mcd1D.scaled(x, h)

  # TODO: Implement a better singularity condition
  if (mcd$var <= .Machine$double.eps^2) {
    candidates = matrix(c(0.0, xbar, xbar, sqrts2, sqrts2, 1.0), nrow = 1L)
  } else {
    # TODO: Let the user choose the threshold
    bulk = x[union(mcd$best, which(abs(x - mcd$mean)/sqrt(mcd$var) <= 2.0))]

    # TODO: Implement a better singularity condition
    if (sd(bulk) <= .Machine$double.eps) {
      candidates = matrix(c(0.0, xbar, xbar, sqrts2, sqrts2, 1.0), nrow = 1L)
    } else {
      left.tail  = x[which(x < min(bulk))]
      right.tail = x[which(x > max(bulk))]

      candidates = matrix(0.0, nrow = 3, ncol = 6)

      # candidate a)
      theta = init.theta(union(left.tail, bulk), right.tail)
      candidates[1, ] = c(theta, 1.0)

      # candidate b)
      theta = init.theta(union(bulk, right.tail), left.tail)
      candidates[2, ] = c(theta, 1.0)

      # candidate c)
      theta = init.theta(bulk, union(left.tail, right.tail))
      candidates[3, ] = c(theta, 1.0)
    }
  }

  candidates = unique(candidates, MARGIN = 1)

  if (nrow(candidates) > 1L) {
    for (i in 1:nrow(candidates)) {
      ML = ML_norm_mix(x, h = h, candidates[i, 1], candidates[i, 2:3], candidates[i, 4:5])
      candidates[i, ] = c(ML$pi, ML$mu, ML$sigma, ML$ward)
    }
  }

  ind = which.min(candidates[, 6])

  return(list(mu.bulk      = if (candidates[ind, 1] <= 0.5) candidates[ind, 2] else candidates[ind, 3],
              sigma.bulk   = if (candidates[ind, 1] <= 0.5) candidates[ind, 4] else candidates[ind, 5],
              ward.measure = candidates[ind, 6]))
}

ML_norm_mix <- function(x, h, pi0, mu0, sigma0) {
  n = length(x)
  b = sqrt(mcd1D.scaled(x, h)$var/n) # artificial lower bound for sigmas

  # negative log-likelihood
  objective  <- function(par) ML_objective(par, x, nargout = 1)$obj
  newton_dir <- function(par) ML_objective(par, x, nargout = 4)$newton.dir
  projector  <- function(par) ML_projector(par, b)

  par0 = c(pi0, mu0, sigma0)

  ftol = 1.0E-6*abs(objective(par0))
  gtol = 1.0E-6*max(abs(newton_dir(par0)))

  #opt = .spg(par = par0, fn = objective, gr = newton_dir, project = projector, method = 3L,
  #                  control = list(maximize = FALSE, maxit = 500, gtol = gtol, ftol = ftol, M = 1), quiet = TRUE, alertConvergence = FALSE)
  opt <- .pgd(par = par0, fn = objective, gr = newton_dir, project = projector,
              control = list(maximize = FALSE, maxit = 500, gtol = gtol, ftol = ftol),
              quiet = TRUE, alertConvergence = FALSE)

  par = opt$par

  pi    = par[1]
  mu    = par[2:3]
  sigma = par[4:5]

  ward = ((1.0 - pi)*sigma[1]^2 + pi*sigma[2]^2)/((1.0 - pi)*sigma[1]^2 + pi*sigma[2]^2 + pi*(1.0 - pi)*(mu[1] - mu[2])^2)

  return(list(pi = pi, mu = mu, sigma = sigma, ward = ward, iter = opt$iter))
}

## R-version:
##
# ML_norm_mix <- function(x, h, pi0, mu0, sigma0) {
#   n = length(x)
#   b = sqrt(mcd1D.scaled(x, h)$var/n) # artificial lower bound for sigmas
#
#   phi       <- function(x, mu, sigma) dnorm(x, mu, sigma, log = FALSE)
#   phi_mu    <- function(x, mu, sigma) ((x - mu)/sigma^2)*phi(x, mu, sigma)
#   phi_sigma <- function(x, mu, sigma) ((x - mu)^2/sigma^2 - 1)/sigma*phi(x, mu, sigma)
#
#   objective <- function(par) {
#     pi = par[1]; mu = par[2:3]; sigma = par[4:5]
#     return(mean(log((1 - pi)*phi(x, mu[1], sigma[1]) + pi*phi(x, mu[2], sigma[2]))))
#   }
#
#   gradient_fisher_inf_matrix <- function(par) {
#     pi       = par[1]; mu = par[2:3]; sigma = par[4:5]
#     f        = ((1 - pi)*phi(x, mu[1], sigma[1]) + pi*phi(x, mu[2], sigma[2]))
#
#     grad_log_f = cbind((-phi(x, mu[1], sigma[1]) + phi(x, mu[2], sigma[2]))/f,
#                        (1 - pi)*phi_mu(x, mu[1], sigma[1])/f,
#                        pi*phi_mu(x, mu[2], sigma[2])/f,
#                        (1 - pi)*phi_sigma(x, mu[1], sigma[1])/f,
#                        pi*phi_sigma(x, mu[2], sigma[2])/f)
#
#     # Empirical Fisher information matrix
#     finf_log_f = crossprod(grad_log_f)/n
#
#     # Gradient
#     grad_log_f = matrix(colMeans(grad_log_f), ncol = 1L)
#
#     return(list(gr = grad_log_f, finf = finf_log_f))
#   }
#
#   projector <- function(par) {
#     pi = par[1]; mu = par[2:3]; sigma = par[4:5]
#     return(c(max(0, min(pi, 1)),
#              mu,
#              pmax(sigma, b)))
#   }
#
#   par0 = c(pi0, mu0, sigma0)
#
#   ftol = 1.0E-6*abs(objective(par0))
#   gtol = 1.0E-6*max(abs(gradient_fisher_inf_matrix(par0)$gr))
#
#   opt <- .fscoring(par = par0, fn = objective, gr_finf = gradient_fisher_inf_matrix, project = projector, method = 3L,
#                    control = list(maximize = TRUE, maxit = 500, gtol = gtol, ftol = ftol, M = 1), quiet = TRUE, alertConvergence = FALSE)
#
#   pi = opt$par[1]
#   mu = opt$par[2:3]
#   sigma = opt$par[4:5]
#
#   # ((1.0 - pi)*sigma1^2 + pi*sigma2^2)/((1.0 - pi)*sigma1^2 + pi*sigma2^2 + pi*(1.0 - pi)*(mu1 - mu2)^2)
#   ward = ((1.0 - pi)*sigma[1]^2 + pi*sigma[2]^2)/((1.0 - pi)*sigma[1]^2 + pi*sigma[2]^2 + pi*(1.0 - pi)*(mu[1] - mu[2])^2)
#
#   return(list(pi = pi, mu = mu, sigma = sigma, ward = ward, iter = opt$iter))
# }
