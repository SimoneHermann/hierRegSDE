#' Bayesian estimation in mixed stochastic differential equations
#'
#' @description Bayesian estimation of the random effects \eqn{\phi_j} in the mixed SDE
#'  \eqn{dY_i(t)= b(\phi_i, t, Y_i(t))dt + \gamma s(t, Y_i(t)) dW_i(t), , \phi_i~N(\mu, \Omega), i=1,...,n} and the parameters
#'  \eqn{\mu, \Omega, \gamma^2}.
#' @param t vector of observation times
#' @param y matrix or list of the n trajectories
#' @param prior list of prior parameters - list(m, v, priorOmega, alpha, beta), priorOmega=list(alpha, beta) if Omega="diag", otherwise prior matrix of Wishart distribution
#' @param start list of starting values
#' @param bSDE b(phi, t, x) drift function
#' @param sVar variance function s^2
#' @param ipred which of the n trajectories is the one to be predicted
#' @param cut the index how many of the ipred-th series are used for estimation
#' @param len number of iterations of the MCMC algorithm - chain length
#' @param Omega structure of the variance matrix Omega of the random effects, diagonal matrix, otherwise inverse wishart distributed
#' @param mod model out of {Gompertz, Richards, logistic, Weibull, Paris, Paris2}, only used instead of bSDE
#' @param modVar default value is sVar(t,x)=1, if "AR": sVar(t,x)=x
#' @param propPar proposal standard deviation of phi is |start$mu|*propPar
#'
#' @return
#' \item{phi}{samples from posterior of \eqn{\phi}}
#' \item{mu}{samples from posterior of \eqn{\mu}}
#' \item{Omega}{samples from posterior of \eqn{\Omega}}
#' \item{gamma2}{samples from posterior of \eqn{\gamma^2}}
#'
#' @examples
#' mod <- "Gompertz"
#' bSDE <- getFun("SDE", mod)
#' mu <- getPar("SDE", mod, "truePar")
#' n <- 5
#' parameters <- defaultPar(mu, n)
#' y <- drawData("SDE", bSDE, parameters)
#' t <- parameters$t
#'
#' prior <- getPrior(mu, parameters$gamma2)
#' start <- getStart(mu, n)
#' chains <- estSDE(t, y, prior, start, bSDE)
#'
#' @details
#' Simulation from the posterior distribution of the random effect from n independent trajectories of the SDE (the Brownian motions \eqn{W1,...,Wn} are independent).
#' \subsection{?}{
#' }
#'
#' @keywords estimation
#' @references Hermann et al. (2015)
#' @export

estSDE <- function(t, y, prior, start, bSDE, sVar, ipred=1, cut, len=1000, Omega="diag", mod=c("Gompertz","logistic","Weibull","Richards","Paris","Paris2"), modVar="", propPar=0.02){
  mod <- match.arg(mod)
  if(is.matrix(y)){
    if(nrow(y)==length(t)){
      y <- t(y)
    }else{
      if(ncol(y)!=length(t)){
        print("length of t has to be equal to the columns of y")
        break
      }
    }
    if(missing(cut)) cut <- length(t)
    t1 <- t
    y1 <- y
    t <- list()
    y <- list()
    for(i in (1:nrow(y1))[-ipred]){
      t[[i]] <- t1
      y[[i]] <- y1[i,]
    }
    t[[ipred]] <- t1[1:cut]
    y[[ipred]] <- y1[ipred,1:cut]
  }
  if(missing(bSDE)) bSDE <- getFun("SDE", mod)
  if(missing(sVar)) sVar <- getFunVar(modVar)

  if(Omega=="diag"){
    postOm <- function(phi,mu){
      postOmegaDiag(prior$priorOmega$alpha,prior$priorOmega$beta,phi,mu)
    }
    if(!is.list(prior$priorOmega)){print("prior parameter for Omega has to be list of alpha and beta")}
  }else{
    postOm <- function(phi,mu){
      postOmega(prior$priorOmega,phi,mu)
    }
    if(!is.matrix(prior$priorOmega)){print("prior parameter for Omega has to be matrix R")}
  }
  propSd <- abs(start$mu)*propPar
  postPhii <- function(lastPhi, mu, Omega, gamma2, X, t){  # X, t vektoren
    lt <- length(t)
    dt <- t[-1]-t[-lt]
    phi_old <- lastPhi
    phi_drawn <- phi_old + rnorm(length(mu),0,propSd)
    ratio <- dmvnorm(phi_drawn,mu,Omega)/dmvnorm(phi_old,mu,Omega)
    ratio <- ratio* prod( dnorm(X[-1], X[-lt] + bSDE(phi_drawn,t[-lt],X[-lt])*dt, sqrt(gamma2*sVar(t[-lt],X[-lt])^2*dt))/dnorm(X[-1], X[-lt] + bSDE(phi_old,t[-lt],X[-lt])*dt, sqrt(gamma2*sVar(t[-lt],X[-lt])^2*dt)))
    if(is.na(ratio)) ratio <- 0
    if(runif(1) <= ratio){
      phi_old <- phi_drawn
    }
    phi_old
  }
  n <- length(y)

  postGamma2 <- function(phi){
    alphaPost <- prior$alpha + sum(sapply(y,length)-1)/2
    help <- numeric(n)
    for(i in 1:n){
      ni <- length(t[[i]])
      delta <- diff(t[[i]])
      help[i] <- sum( (y[[i]][-1] - y[[i]][-ni] - bSDE(phi[i,],t[[i]][-ni],y[[i]][-ni])*delta)^2/(sVar(t[[i]][-ni],y[[i]][-ni])^2*delta) )
    }
    betaPost <-  prior$beta + sum(help)/2
    1/rgamma(1, alphaPost, betaPost)
  }

  phi_out <- list()
  mu_out <- matrix(0,len,length(start$mu))
  Omega_out <- list()
  gamma2_out <- numeric(len)

  phi <- start$phi
  gamma2 <- start$gamma2
  mu <- start$mu
  Omega <- postOm(phi,mu)

  for(count in 1:len){

    for(i in 1:n){
      phi[i,] <- postPhii(phi[i,], mu, Omega, gamma2, y[[i]], t[[i]])
    }
    mu <- postmu(phi, prior$m, prior$v, Omega)
    Omega <- postOm(phi, mu)
    gamma2 <- postGamma2(phi)

    phi_out[[count]] <- phi
    mu_out[count,] <- mu
    Omega_out[[count]] <- Omega
    gamma2_out[count] <- gamma2
  }
  list(phi=phi_out, mu=mu_out, Omega=Omega_out, gamma2=gamma2_out)
}
