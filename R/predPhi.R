#' Prediction Of The Random Effects In Mixed Stochastic Differential Equations
#'
#' @description Prediction of the random effects
#' @param samples output of estSDE or estREg
#' @param cand candidates for phi (matrix with p columns)
#' @return matrix phi
#'
#'

#' @examples
#'
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
#' phi.pred <- predPhi(chains)
#'
#'
#' @keywords prediction
#' @references
#' follows...
#'
#' @export
predPhi <- function(samples, cand){ # only for the case Omega="diag"

  K <- nrow(samples$mu)
  p <- ncol(samples$mu)
  if(missing(cand)){
    cand <- sapply(1:p, function(i){
    seq(mean(samples$mu[,i])-4*sqrt(mean(sapply(samples$Omega, function(mat) mat[i,i]))), mean(samples$mu[,i])+4*sqrt(mean(sapply(samples$Omega, function(mat) mat[i,i]))), length=500)
  })
  }
  phi.res <- matrix(0, K, p)
  for(i in 1:p){
    prob <- numeric(length(cand[,i]))
    for(a in 1:length(cand[,i])){
      prob[a] <- mean(dnorm(cand[a,i], samples$mu[,i], sqrt(sapply(samples$Omega, function(mat) mat[i,i]))))
    }
    cand2 <- cand[prob!=0,i]
    prob <- prob[prob!=0]
    mp <- max(prob)
    count <- 1
    while(count <= K){
      u <- runif(1,0,mp)
      ind <- sample(1:length(cand2),1)
      if(u <= prob[ind]){
        phi.res[count,i] <- cand2[ind]
        count <- count + 1
      }
    }
  }
  return(phi.res)
}

