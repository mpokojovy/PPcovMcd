## MCD with asympotically unbiased scaling
mcd1D.scaled <- function(x, h) {
  mcd  = mcd1D(x, h)
  frac = h/length(x)
  c    = frac/pchisq(qchisq(frac, df = 1), df = 3)
  return(list(mean = mcd$mean, var = c*mcd$var, best = mcd$best, h = h))
}

## Matrix multiplication with arma library
"%arma*%" <- function(A, B) arma_ast(A, B)

## Mahalanobis standartization
mah.standard <- function(x, loc = NULL, cov = NULL, rank = ncol(x)) {
  if (is.null(loc)) loc = base::colMeans(x)
  if (is.null(cov)) cov = stats::cov(x)

  return(mah_standard(x, loc, cov, rank))
}