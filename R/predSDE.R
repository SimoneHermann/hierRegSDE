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
#' chains <- estSDE(t, y, prior, cut = 50, ipred = 1, start, bSDE, len = 5000)
#' ind <- seq(1001, 5000, by = 4)
#' samples <- list(phi = sapply(1:length(mu), function(i) phi = phi_ij(chains$phi, 1, i))[ind, ], gamma2 = chains$gamma2[ind])
#' prediction <- predSDE(t[50:101], samples, y[1,50], bSDE, cand = seq(-2, 2, length = 1000))
#' plot(t[51:101], y[1,51:101], ylim = range(y[1,51:101]) + c(0, 1))
#' lines(t[51:101], apply(prediction, 1, quantile, 0.025), col = 3)
#' lines(t[51:101], apply(prediction, 1, quantile, 0.975), col = 3)
#' lines(t[51:101], apply(prediction, 1, mean), col = 2)
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
