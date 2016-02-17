#' Bayesian prediction in mixed SDE models
#'
#' @description Bayesian prediction in the mixed SDE
#'  \eqn{dY_j(t)= b(\phi_j, t, Y_j(t))dt + \gamma s(t, Y_j(t)) dW_j(t), , \phi_j~N(\mu, \Omega)}.
#' @param tau vector of times which are predicted
#' @param samples list of samples from the posterior
#' @param last last observation - staring point for the prediction of Markov chain
#' @param bSDE drift function
#' @param sVar variance function
#' @param cand vector of candidates for trajection sampling
#' @param mod model out of {Gompertz, Richards, logistic, Weibull}, only used instead of bSDE
#' @param modVar model for the variance structure
#'
#' @return matrix of predictions in t
#' @export



predSDE <- function(tau, samples, last, bSDE, sVar, cand, mod="Gompertz", modVar=""){
  # tau = (tni,t^*), where t^* can be a vector
  if(missing(bSDE)) bSDE <- getFun("SDE", mod)
  if(missing(sVar)) sVar <- getFunVar(modVar)
  if(missing(cand)) cand <- seq(-1,1,by=0.01)

  K <- length(samples$gamma2)
  if(ncol(samples$phi)==K) samples$phi <- t(samples$phi)

  L <- length(tau)-1
  Dt <- tau[-1]-tau[-(L+1)]
  Xtau <- matrix(0,L,K)
  vec <- rep(last,K)
  prob <- matrix(0,K,length(cand))
  cand_l <- cand + last

  for(l in 1:L){
    for(a in 1:K){
      prob[a,] <- dnorm(cand_l, vec[a]+bSDE(samples$phi[a,],tau[l],vec[a])*Dt[l], sqrt(samples$gamma2[a]*Dt[l])*sVar(tau[l],vec[a]))
    }

    prEnd <- apply(prob,2,mean)
    cand_l <- cand_l[prEnd!=0]
    prEnd <- prEnd[prEnd!=0]

    mp <- max(prEnd)

    count <- 1
    while(count <= K){
      u <- runif(1,0,mp)
      ind <- sample(1:length(cand_l),1)
      if(u <= prEnd[ind]){
        vec[count] <- cand_l[ind]
        count <- count + 1
      }
    }
    cand_l <- cand + mean(vec)
    Xtau[l,] <- vec
  }
  Xtau
}
