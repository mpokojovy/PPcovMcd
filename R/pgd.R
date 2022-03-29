## BB::spg rewritten as Projected Gradient Descent with Armijo's Rule
## https://rdrr.io/rforge/BB/src/R/spg.R

.pgd <- function(par, fn, gr, project, control = list(), quiet = FALSE, alertConvergence = TRUE) {
  # control defaults
  ctrl <- list(maximize = TRUE, maxit = 1500, ftol = 1.0E-10, gtol = 1.0E-05, maxfeval = 10000L, trace = !quiet, triter = 10L)
  ctrl[intersect(names(ctrl), names(control))] <- control[intersect(names(ctrl), names(control))]
  
  maxit    <- ctrl$maxit
  ftol     <- ctrl$ftol
  gtol     <- ctrl$gtol
  maxfeval <- ctrl$maxfeval
  maximize <- ctrl$maximize
  trace    <- ctrl$trace
  triter   <- ctrl$triter
  
  func <- if (maximize) function(par) -fn(par) else fn
  grad <- if (maximize) function(par) -gr(par) else gr
  
  gamma = 1E-4
  
  p = length(par)
  
  armijo <- function(par, fv, step, dir, feval, func) {
    lasteq = 0L
    
    lastpar  = par
    lastfv   = fv
    laststep = step
    newstep  = step
    
    repeat {
      newpar = project(par + newstep*dir)
      
      newfv  = func(newpar)
      feval  = feval + 1L
      
      if (feval > maxfeval)
        return(list(par = NA, f = NA, feval = NA, lsflag = 2))
      
      dpar = newpar - par
      
      neweq = if ((newfv - fv) <= -gamma/newstep*sum(dpar*dpar)/p) 1L else -1L
      
      if ((lasteq >= 0) && (neweq > 0)) {
        if (newstep >= 1.0) {
          break
        } else {
          lastpar  = newpar
          lasteq   = neweq
          lastfv   = newfv
          laststep = newstep
          
          newstep = laststep*2.0
        }
      } else if ((lasteq <= 0) && (neweq < 0)) {
        lastpar  = newpar
        lasteq   = neweq
        lastfv   = newfv
        laststep = newstep
        
        newstep = laststep*0.5
      } else if ((lasteq < 0) && (neweq > 0)) {
        break 
      } else if ((lasteq > 0) && (neweq < 0)) {
        newpar  = lastpar
        newfv   = lastfv
        newstep = laststep
        break
      }
    }
    
    return(list(par = newpar, f = newfv, step = newstep, feval = feval, lsflag = 0L))
  }
  
  step = 1.0
  
  par = project(par)
  g   = grad(par)
  dir = step*g
  
  f0    = f = func(par)
  feval = 1L
  
  fchg    = Inf
  parchg2 = Inf
  
  iter = 0L
  
  if (trace) cat("iter: ", 0L, " f-value: ", f0*(-1)^maximize, " step: ", 1.0, "\n")
  
  while ((iter <= maxit) && (fchg > ftol) && (parchg2 > gtol)) {
    if (is.infinite(sum(g*g/p))) {
      lsflag = 4L
      break
    }
    
    iter = iter + 1L
    
    armijo.out = armijo(par = par, fv = f, step = step, dir = -g, feval = feval, func = func)
    
    lsflag = armijo.out$lsflag
    
    if (lsflag != 0) break
    
    fchg   = abs(f - armijo.out$f)
    f      = armijo.out$f
    step   = armijo.out$step
    feval  = armijo.out$feval
    
    parnew = armijo.out$par
    gnew   = grad(parnew)
    
    parchg  = parnew - par
    parchg2 = sqrt(sum(parchg*parchg)/p)
    
    par = parnew
    g   = gnew
    
    if (trace && (iter %% triter == 0))
      cat("iter: ", iter, " f-value: ", f*(-1)^maximize, " step: ", step, "\n")
  }
  
  if (is.null(lsflag)) {
    if (!quiet) warning("convergence tolerance satisified at intial parameter values.")
    lsflag <- 0
  }
  
  type = 0 # "Successful convergence")
  if (lsflag == 0) {
    #if ((fchg <= ftol) && (parchg2 <= gtol)) type = 0 # "Successful convergence")
    if (iter >= maxit) type = 1 # "Maximum number of iterations exceeded")
  } else {
    if (lsflag == 1) type = 3 # "Failure:  Error in function evaluation")
    if (lsflag == 2) type = 2 # "Maximum function evals exceeded")
    if (lsflag == 3) type = 4 # "Failure:  Error in gradient evaluation")
    if (lsflag == 4) type = 5 # "Failure:  Error in projection")
  }
  
  if (alertConvergence && (type != 0))
    warning("Unsuccessful convergence.")
  
  return(list(par = par, value = (-1)^maximize*f, gradient = (-1)^maximize*g, 
              fn.reduction = (-1)^maximize*(f0 - f), 
              iter = iter, feval = feval, convergence = type))
}  