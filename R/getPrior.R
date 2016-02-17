#' Builds list of prior parameters
#'
#' @description Creation of list of parameters conditional on true values \eqn{\mu} and \eqn{\gamma^2}
#' @param mu
#' @param gamma2
#' @param mixed
#' @return list of prior values
#' @export
getPrior <- function(mu, gamma2, mixed=TRUE){
  if(mixed){
    prior <- list( m=mu, v=rep(1,length(mu)),priorOmega=list(alpha=rep(100,length(mu)), beta=rep(1,length(mu))), alpha=1/gamma2+1, beta=1)
  }else{
    prior <- list(mu=mu, Omega=mu, alpha=1/gamma2+1, beta=1)
  }
  prior
}
