#' Posterior
#'
#' @description Posterior for parameters \eqn{\Omega}
#' @param alpha
#' @param beta
#' @param phi
#' @param mu
#' @return one sample of posterior


postOmegaDiag <- function(alpha, beta, phi, mu){  # length(alpha)=length(beta)=length(mu)
  p <- length(mu)
  Dia <- numeric(p)
  for(i in 1:p){
    Dia[i] <- 1/rgamma(1, alpha[i] + nrow(phi)/2, beta[i] + sum((phi[,i]-mu[i])^2)/2)
  }
  diag(Dia)
}
