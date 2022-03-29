## M34 directions
.M34_step <- function(x, xbar = colMeans(x), S = cov(x), h = floor((nrow(x) + ncol(x) + 1)/2)) {
  n = nrow(x); p = ncol(x)

  # 3 candidates (eigen)directions as column vectors
  v = matrix(0.0, nrow = p, ncol = 3)

  M34 = M34_vec(x, xbar, S)

  v[, 1] = M34$v1
  v[, 2] = M34$v2
  v[, 3] = M34$v3
  rm(M34)

  mu    = c(0.0, 0.0, 0.0)
  sigma = c(0.0, 0.0, 0.0)
  ward  = c(0.0, 0.0, 0.0)

  var.ratio = c(0.0, 0.0, 0.0)

  for (i in 1:3) {
    xp = x %*% v[, i] # projection
    var.ratio[i] = mcd1D(xp, h)$var
  }

  var.ratio = var.ratio/max(var.ratio)

  for (i in 1:3) {
    xp = x %*% v[, i] # projection

    # TODO: Implement a better smallness test
    if (var.ratio[i] < .Machine$double.eps) {
      if (i != 3) {
        # The second candidate direction is degerate for the projected data
        mu[i]    = NA
        sigma[i] = NA
        ward[i]  = NA
      } else {
        mu[i]    = mean(xp)
        sigma[i] = sd(xp)
        ward[i]  = 1.0
      }
    } else {
      ML = ML_bulk(xp, h)

      mu[i]    = ML$mu.bulk
      sigma[i] = ML$sigma.bulk
      ward[i]  = ML$ward.measure
    }
  }

  ind = max(which(abs(ward - min(ward, na.rm = TRUE)) < sqrt(.Machine$double.eps)))

  # ward.measure = bimodality projection index
  return(list(dir = as.matrix(v[, ind], ncol = 1L),
              mu.bulk = mu[ind],
              sigma.bulk = sigma[ind],
              ward.measure = ward[ind]))
}

## Apply gradient descend to converge at robust principal direction starting at initial guess
.PP_step <- function(x, init.dir, h = floor((nrow(x) + ncol(x) + 1)/2), trace = as.integer(trace), maxiter = 100L) {
  n = nrow(x); p = ncol(x)

  dir.old  = init.dir
  best.old = as.vector(mcd1D(x %*% init.dir, h)$best)

  iter = 0L

  while (iter <= maxiter) {
    objective <- function(dir) PP_objective(x %*% dir, best.old)
    gradient  <- function(dir) PP_gradient(x %*% dir, x, best.old)
    projector <- function(dir) PP_project(dir)

    ftol = abs(objective(dir.old))*1E-6
    gtol = sqrt(sum(gradient(dir.old)^2/p))*1E-6

    opt <- .pgd(par = dir.old, fn = objective, gr = gradient, project = projector,
                control = list(maximize = TRUE, gtol = gtol, ftol = ftol),
                quiet = (trace < 1L), alertConvergence = FALSE)

    dir.new  = opt$par
    best.new = as.vector(mcd1D(x %*% dir.new, h)$best)

    iter = iter + 1L

    if (setequal(best.old, best.new)) {
      break
    }

    dir.chg = dir.new - dir.old

    if (sqrt(sum(dir.chg*dir.chg)/p) <= 1E-9) break

    dir.old  = dir.new
    best.old = best.new
  }

  return(dir.new)
}

##' Implements Pokojovy and Jobe's (2020) robust projection pursuit (PP) method to obtain a robust initial set
##'
##' @title Robust initial observation ordering based on Jobe and Pokojovy's (2019) robust projection-pursuit method
##' @param x  $(n \times p)$-dimensional data matrix.
##' @param h  bulk size integer. Any value between $\floor{n + p + 1}{2}\rfloor$ and $n$ are allowed.
##' @param full.h Return full (length $n$) ordering or only the first $h$ observation?
##' @param trace Display gradient descent progress?
##' @return Optimal direction as a column vector of length $k$
##' @author Michael Pokojovy
.pp_initset <- function(x, h = floor((nrow(x) + ncol(x) + 1)/2), full.h = FALSE, trace = as.integer(trace), standardization = "none") {
  n = nrow(x); p = ncol(x)

  if (p == 1L) {
    if (trace > 0L)
      cat("Applying exact 1D MCD algorithm...\n")

    return(mcd1D(x = x, h = h)$best)
  }

  if (standardization == "spherical") {
    # Apply spherical standardization
    x = mah.standard(x)$z.scores
  } else if (standardization == "diagonal") {
    ## Apply diagonal standardization
    for (i in 1:p) {
      Qn.hat = robustbase::s_Qn(x[, i], mu.too = TRUE)
      x[, i] = (x[, i] - Qn.hat[1])/Qn.hat[2]
    }
  }

  ## Determine principal directions
  prin.dir = matrix(0.0, nrow = p, ncol = p)
  xp = x # iteratively projected/reduced dataset
  locs = rep(0.0, p) # robust location
  vars = rep(0.0, p) # robust principal variances

  norm.dim.cnt  = 0L    # number of sequential projections looking (more or less) normal

  dim = 0L
  repeat {
    dim = dim + 1L
    if (dim >= p) {
      break
    }

    if (trace > 0L) {
      cat("PP method: principal direction # =", dim, " ... computing \n")
    }

    ## Obtain next robust principal direction
    if ((dim < p) || (nrow(xp) > h)) {
      # Guess initial direction
      M34.est = .M34_step(xp, xbar = colMeans(xp),
                          S = if (dim == 1) cov(xp) else cov(xp) + tcrossprod(prin.dir[, 1:(dim - 1)]), h)

      dir = M34.est$dir

      dir1 = dir
      ward.measure1 = M34.est$ward.measure
      rm(M34.est)

      # Try to improve initial guess
      dir = .PP_step(xp, init.dir = dir, h = h, trace = trace)

      dir2 = dir
      ward.measure2 = ML_bulk(xp %*% dir2, h)$ward.measure

      dir = if ((ward.measure1 - ward.measure2) <= sqrt(.Machine$double.eps)) dir1 else dir2

      y = as.vector(xp %*% dir)
      mcd = mcd1D.scaled(y, h)

      if (sqrt(mcd$var) <= .Machine$double.eps) {
        break
      }

      locs[dim] = mcd$mean
      vars[dim] = mcd$var

      #TODO: Let user choose thresholds
      norm.dim.cnt = if (var(y)/vars[dim] >= 1.25) 0L else norm.dim.cnt + 1L
      if (norm.dim.cnt >= max(0.05*p, 5L)) {
        prin.dir[, dim] = dir
        dim = dim + 1L
        break
      }
    }
    else {
      break
    }

    # TODO: Improve singularity test
    if (sqrt(vars[dim]) <= .Machine$double.eps) {
      stop("Data are numerically not of full rank.")
    }

    prin.dir[, dim] = dir

    z = abs(y - locs[dim])/sqrt(vars[dim]) # compute robust z-scores

    # TODO: Let user select respective normal quantile
    bulk = union(order(z)[1:h], which(z <= 2.0))

    xp = xp[bulk, ]

    if (length(bulk) == h) {
      dim = dim + 1L
      break
    }

    ## Project the data onto the orthogonal space
    xp = proj_orth_comp(xp, dir)
  }

  ## Add missing directions
  if (dim < p) {
    eig = eigen(M4_matrix(xp, xbar = colMeans(xp),
                          S = if (dim == 1) cov(xp) else cov(xp) + tcrossprod(prin.dir[, 1:(dim - 1)])))
    prin.dir[, dim:p] = eig$vectors[, 1:(p - dim + 1)]
    rm(eig)
  } else {
    prin.dir[, p] = pracma::nullspace(t(prin.dir[, 1:(p - 1)]))[, 1L]
  }

  if (dim <= p) {
    for (d in dim:p) {
      if (trace > 0L) {
        cat("PP method: principal direction # =", d, " ... computing \n")
      }

      y = as.vector(xp %*% prin.dir[, d])

      if (nrow(xp) > h) {
        mcd = mcd1D.scaled(y, h)

        locs[d] = mcd$mean
        vars[d] = mcd$var
      } else {
        locs[d] = mean(y)
        vars[d] = var(y)
      }

      # TODO: Improve singularity test
      if (sqrt(vars[d]) <= .Machine$double.eps) {
        stop("Data are numerically not of full rank.")
      }
    }
  }

  ## Compute robust squared Mahalanobis distances
  xp = x %*% prin.dir # since prin.dir is "orthogonal", this corresponds to linear transformation v = S^{-1} u, u column vector

  for (dim in 1L:p) {
    xp[, dim] = (xp[, dim] - locs[dim])/sqrt(vars[dim]) # standardize the i-th robust principal component
  }

  d2 = rowSums(xp^2)

  return(if (full.h) sort.list(d2) else sort.list(d2)[1L:h])
}
