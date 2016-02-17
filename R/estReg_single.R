#' Bayesian stimation in nonlinear regression models
#'
#' @description Bayesian estimation of the parameters of the mixed nonlinear regression model
#'  \eqn{y_j= f(\phi, t_j) + \epsilon_j, \epsilon_j~N(0,\gamma^2*s^2(t_j)}.
#' @param t vector of observation times
#' @param y vector of the M trajectories
#' @param prior list of prior parameters - list(mu, Omega, alpha, beta)
#' @param start list of starting values
#' @param fODE regression function
#' @param sVar variance function
#' @param len number of iterations of the MCMC algorithm
#' @param mod model out of {Gompertz, Richards, logistic, Weibull}, only used instead of fODE
#' @param propSd proposal standard deviation of phi is |mu|*propSd
#'
#' @return
#' \item{phi}{estimator of \eqn{\phi}}
#' \item{gamma2}{estimator of \eqn{\gamma^2}}
#' @export

estReg_single <- function(t, y, prior, start, fODE, sVar, len=1000, mod="Gompertz"){  # prior liste von prior-parametern, len=Anzahl von Ziehungen

  if(missing(fODE)) fODE <- getFun("ODE", mod)
  if(missing(sVar)) sVar <- getFunVar()

  propSd <- abs(prior$mu)/50
  lt <- length(t)
  lphi <- length(start$phi)

  #   postPhi <- function(lastPhi, gamma2){
  #     phi_old <- lastPhi; phi_drawn <- lastPhi
  #     for(k in 1:lphi){
  #       phi_drawn[k] <- rnorm(1, phi_old[k], propSd[k])
  #       ratio <- prod(dnorm(y, fODE(phi_drawn,t), sqrt(gamma2*sVar(t)))/dnorm(y, fODE(phi_old,t), sqrt(gamma2*sVar(t))))
  #       ratio <- ratio*dnorm(phi_drawn[k], prior$mu[k], prior$Omega[k])
  #       if(is.na(ratio)) ratio <- 0
  #       if(runif(1)<ratio){
  #         phi_old[k] <- phi_drawn[k]
  #       }else{
  #         phi_drawn[k] <- phi_old[k]
  #         }
  #     }
  #     phi_old
  #   }
  postPhi <- function(lastPhi, gamma2){
    phi_old <- lastPhi
    phi_drawn <- rnorm(lphi,phi_old,propSd)
    ratio <- dmvnorm(phi_drawn,prior$mu,diag(prior$Omega))/dmvnorm(phi_old,prior$mu,diag(prior$Omega))
    ratio <- ratio* prod(dnorm(y, fODE(phi_drawn,t), sqrt(gamma2*sVar(t)))/dnorm(y, fODE(phi_old,t), sqrt(gamma2*sVar(t))))
    if(is.na(ratio)){ratio <- 0}
    if(runif(1)<ratio){
      phi_old <- phi_drawn
    }
    phi_old
  }

  postGamma2 <- function(phi){
    alphaPost <- prior$alpha + lt/2
    betaPost <-  prior$beta + sum((y-fODE(phi, t))^2/sVar(t))/2
    1/rgamma(1, alphaPost, betaPost)
  }

  phi_out <- matrix(0,len,length(prior$mu))
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
