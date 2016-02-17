#' Posterior
#'
#' @description Posterior for parameters \eqn{\Omega}
#' @param R
#' @param phi
#' @param mu
#' @return one sample of posterior



postOmega <- function(R, phi, mu){
  Rpost <- solve(R + (t(phi)-as.vector(mu))%*%t((t(phi)-as.vector(mu))))
  solve( rWishart(1,nrow(phi)+length(mu)+1,Rpost)[,,1])
}
