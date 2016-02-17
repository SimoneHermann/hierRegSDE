#' Bayesian stimation in mixed stochastic differential equations
#'
#' @description Bayesian estimation of the parameters in the mixed SDE
#'  \eqn{dY(t)= b(\phi, t, Y(t))dt + \gamma s(t, Y(t)) dW(t)}.
#' @param t vector of observation times
#' @param X vector of the M trajectories
#' @param prior list of prior parameters - list(mu, Omega, alpha, beta)
#' @param start list of starting values
#' @param fODE regression function
#' @param sVar variance function
#' @param len number of iterations of the MCMC algorithm
#' @param mod model out of {Gompertz, Richards, logistic, Weibull}, only used instead of fODE
#' @param modVar default value is sVar(t,x)=1, if "AR": sVar(t,x)=x
#' @param propSd proposal standard deviation of phi is |mu|*propSd
#'
#' @return
#' \item{phi}{estimator of \eqn{\phi}}
#' \item{gamma2}{estimator of \eqn{\gamma^2}}
#' @export

estSDE_single <- function(t, X, prior, start, bSDE, sVar, len=1000, mod="Gompertz", modVar=""){  # p liste von prior-parametern, len=Anzahl von Ziehungen

  if(missing(bSDE)) bSDE <- getFun("SDE", mod)
  if(missing(sVar)) sVar <- getFunVar(modVar)

  propSd <- abs(prior$mu)/50
  lt <- length(t)
  dt <- t[-1]-t[-lt]
  lphi <- length(propSd)

  postPhi <- function(lastPhi, gamma2){
    phi_old <- lastPhi

    phi_drawn <- phi_old + rnorm(lphi,0,propSd)
    ratio <- dmvnorm(phi_drawn,prior$mu,diag(prior$Omega))/dmvnorm(phi_old,prior$mu,diag(prior$Omega))
    ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn,t[-lt],X[-lt])*dt, sqrt(gamma2*sVar(t[-lt],X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old,t[-lt],X[-lt])*dt, sqrt(gamma2*sVar(t[-lt],X[-lt])^2*dt)))
    if(is.na(ratio)){ratio <- 0}
    if(runif(1)<ratio){
      phi_old <- phi_drawn
    }
    phi_old
  }
  postGamma2 <- function(phi){
    alphaPost <- prior$alpha + (lt-1)/2
    betaPost <-  prior$beta + sum( (X[-1] - X[-lt] - bSDE(phi,t[-lt],X[-lt])*dt)^2/(sVar(t[-lt],X[-lt])^2*dt) )/2
    1/rgamma(1, alphaPost, betaPost)
  }

  phi_out <- matrix(0,len,lphi)
  gamma2_out <- numeric(len)

  phi <- start$phi
  gamma2 <- start$gamma2

  for(count in 1:len){

    phi <- postPhi(phi, gamma2)
    gamma2 <- postGamma2(phi)

    phi_out[count,] <- phi
    gamma2_out[count] <- gamma2

  }
  list(phi=phi_out, gamma2=gamma2_out)
}
