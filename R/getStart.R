#' Builds list of starting values
#'
#' @description Creation of list of parameters conditional on true value \eqn{\mu}
#' @param mu
#' @param n
#' @param mixed
#' @return list of starting values
#' @export

getStart <- function(mu, n=1, mixed=TRUE){
  if(mixed){
    start <- list( phi=matrix(rep(mu,each=n),n), gamma2=1, mu=mu)
  }else{
    start <- list(phi=mu, gamma2=1)
  }
  start
}
