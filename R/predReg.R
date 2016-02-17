#' Bayesian prediction in mixed nonlinear regression models
#'
#' @description Bayesian prediction in the mixed nonlinear regression model
#'  \eqn{y_{ij}= f(\phi_j, t_{ij}) + \epsilon_{ij}, \epsilon_{ij}~N(0,\gamma^2*s^2(t_{ij}), \phi_j~N(\mu, \Omega)}.
#' @param t vector of times which are predicted
#' @param samples list of samples from the posterior
#' @param fODE regression function
#' @param sVar variance function
#' @param cand vector of candidates for trajection sampling
#' @param len number of samples from the predictive distribution
#' @param mod model out of {Gompertz, Richards, logistic, Weibull}, only used instead of fODE
#'
#' @return matrix of predictions in t
#' @export


predReg <- function(t, samples, fODE, sVar, cand, len=1000, mod="Gompertz"){

  if(missing(fODE)) fODE <- getFun("ODE", mod)
  if(missing(sVar)) sVar <- getFunVar()
  K <- length(samples$gamma2)
  if(ncol(samples$phi)==K) samples$phi <- t(samples$phi)

  if(missing(cand)){
    he <- fODE(apply(samples$phi, 2, mean), mean(t))
    cand <- seq(he/2, he*2, length=100)
  }

  pred <- function(tSt){
    mat <- matrix(0,length(cand), K)
    for(i in 1:K){
      mat[,i] <- dnorm(cand, fODE(samples$phi[i,],tSt), sqrt(samples$gamma2[i]*sVar(tSt)))
    }
    prob <- apply(mat,1,mean)
    cand <- cand[prob>0]
    prob <- prob[prob>0]
    mp <- max(prob)
    vec <- numeric(len)

    count <- 1
    while(count<=len){
      u <- runif(1,0,mp)
      ind <- sample(1:length(cand),1)
      if(u<=prob[ind]){
        vec[count] <- cand[ind]
        count <- count + 1
      }
    }
    vec
  }
  sapply(t, pred)
}
