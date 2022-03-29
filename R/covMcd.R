## Portions of code are adopted from robustbase R package by Valentin Todorov

## This program is a free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## http://www.r-project.org/Licenses/

.getDefCtrl <- function(nm, defCtrl = rrcov.control()) {
  callerEnv <- parent.frame()
  if(is.null(get(nm, envir = callerEnv)))
    assign(nm, defCtrl[[nm]], envir=callerEnv)
}

PPcovMcd <- function(x,
                   cor = FALSE,
                   raw.only = FALSE,
                   alpha = control$ alpha,
                   nsamp = control$ nsamp,
                   nmini = control$ nmini, kmini = control$ kmini,
                   scalefn=control$scalefn, maxcsteps=control$maxcsteps,
                   initHsets = NULL, save.hsets = FALSE, names = TRUE,
                   seed  = control$ seed,
                   tolSolve = control$ tolSolve, # had 1e-10 hardwired {now 1e-14 default}
                   trace = control$ trace,
                   use.correction = control$ use.correction,
                   wgtFUN = control$ wgtFUN,
                   control = robustbase::rrcov.control(),
                   PP.standardization = "spherical")
{
  if (!is.element(tolower(nsamp), c("pp", "pp.plus"))) {
    mcd = robustbase::covMcd(x = x, cor = cor, raw.only = raw.only, alpha = alpha,
                             nsamp = nsamp, nmini = nmini, kmini = kmini,
                             scalefn = scalefn, maxcsteps = maxcsteps,
                             initHsets = initHsets, save.hsets = save.hsets, names = names,
                             seed  = seed, tolSolve = tolSolve,
                             trace = trace, use.correction = use.correction, wgtFUN = wgtFUN, control = control)

    return(mcd)
  } else {
    if (is.null(initHsets)) {
      if (is.data.frame(x))
        x.copy <- data.matrix(x, rownames.force = FALSE)
      else if (!is.matrix(x))
        x.copy <- matrix(x, length(x), 1,
                         dimnames = if(names) list(names(x), deparse(substitute(x))))
      else
        x.copy <- x

      dimnames(x.copy) <- NULL

      ## drop all rows with missing values (!!) :
      ok <- is.finite(x.copy %*% rep.int(1, ncol(x)))
      x.copy <- x.copy[ok, , drop = FALSE]
      if (!length(dx <- dim(x.copy)))
        stop("All observations have missing values!")
      n <- dx[1]; p <- dx[2]

      dx <- dim(x); n <- dx[1]; p <- dx[2]
      stopifnot(p >= 1, n >= 1)
      h <- robustbase::h.alpha.n(alpha, n, p)

      initHsets = .pp_initset(x.copy, h = h, full.h = FALSE, trace = trace, standardization = PP.standardization)
    } else {
      stop("InitSets is expected to be NULL")
    }

    nsamp = "deterministic" # To make sure no internal problems with covMcd() arise
    ppmcd = robustbase::covMcd(x = x, cor = cor, raw.only = raw.only, alpha = alpha,
                               nsamp = nsamp, nmini = nmini, kmini = kmini,
                               scalefn = scalefn, maxcsteps = maxcsteps,
                               initHsets = initHsets,
                               save.hsets = save.hsets, names = names,
                               seed  = seed, tolSolve = tolSolve,
                               trace = trace, use.correction = use.correction, wgtFUN = wgtFUN, control = control)

    if (tolower(nsamp) == "pp") {
      ppmcd$nsamp = "PP"
      ppmcd$method = gsub("Deterministic", "PP", mcd$method)
      ppmcd$iBest = NA

      return(ppmcd)
    } else {
      nsamp = "deterministic" # To make sure no internal problems with covMcd() arise
      detmcd = robustbase::covMcd(x = x, cor = cor, raw.only = raw.only, alpha = alpha,
                                  nsamp = nsamp, nmini = nmini, kmini = kmini,
                                  scalefn = scalefn, maxcsteps = maxcsteps,
                                  initHsets = NULL,
                                  save.hsets = save.hsets, names = names,
                                  seed  = seed, tolSolve = tolSolve,
                                  trace = trace, use.correction = use.correction, wgtFUN = wgtFUN, control = control)

      if (ppmcd$crit <= detmcd$crit) {
        mcd = ppmcd

        if (ppmcd$crit < detmcd$crit)
         mcd$iBest = 7
        else
         mcd$iBest = c(detmcd$iBest, 7)
      } else {
        mcd = detmcd
      }

      mcd$nsamp = "PP.plus"
      mcd$method = gsub("Deterministic", "PP+", mcd$method)

      return(mcd)
    }
  }


  # } else {
  #   if (is.null(initHsets)) {
  #
  #     if (is.data.frame(x))
  #       x.copy <- data.matrix(x, rownames.force = FALSE)
  #     else if (!is.matrix(x))
  #       x.copy <- matrix(x, length(x), 1,
  #                        dimnames = if(names) list(names(x), deparse(substitute(x))))
  #
  #     dimnames(x.copy) <- NULL
  #
  #     ## drop all rows with missing values (!!) :
  #     ok <- is.finite(x.copy %*% rep.int(1, ncol(x)))
  #     x.copy <- x.copy[ok, , drop = FALSE]
  #     if (!length(dx <- dim(x.copy)))
  #       stop("All observations have missing values!")
  #     n <- dx[1]; p <- dx[2]
  #
  #     dx <- dim(x); n <- dx[1]; p <- dx[2]
  #     stopifnot(p >= 1, n >= 1)
  #     h <- robustbase::h.alpha.n(alpha, n, p)
  #
  #     initHsets = .pp_initset(x.copy, h = h, full.h = FALSE, trace = trace, standardization = PP.standardization)
  #   }
  #
  #   nsamp = "deterministic" # To make sure no internal problems with covMcd() arise
  #
  #   mcd = robustbase::covMcd(x = x, cor = cor, raw.only = raw.only, alpha = alpha,
  #                            nsamp = nsamp, nmini = nmini, kmini = kmini,
  #                            scalefn = scalefn, maxcsteps = maxcsteps,
  #                            initHsets = initHsets,
  #                            save.hsets = save.hsets, names = names,
  #                            seed  = seed, tolSolve = tolSolve,
  #                            trace = trace, use.correction = use.correction, wgtFUN = wgtFUN, control = control)
  #
  #   mcd$nsamp = "PP"
  #   mcd$method = gsub("Deterministic", "PP", mcd$method)
  #
  #   return(mcd)
  # }
} ## {covMcd}
