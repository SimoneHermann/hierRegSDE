#' Default parameters
#'
#' @description Yields example parameters conditional on \eqn{\mu} for the different models.
#' @param mu vector of paramters
#' @param n number of series (if mixed=TRUE)
#' @param mixed default=TRUE; if false: \eqn{\phi_i=\mu}
#' @return
#' \item{phi}{example values of \eqn{\phi}}
#' \item{mu}{example values of \eqn{\mu}}
#' \item{Omega}{example values of \eqn{\Omega}}
#' \item{gamma2}{example values of \eqn{\gamma^2}}
#' \item{t}{example time points}
#' \item{startY}{starting value for the SDE}
#' @export

defaultPar <- function(mu, n=20, mixed=TRUE){
  Omega <- numeric(length(mu))
  Omega[mixed] <- abs(mu)/100
  phi <- t(replicate(n, rnorm(length(mu), mu, sqrt(Omega))))
  gamma2 <- (1/2)^2
  t <- seq(0,1,length=101)

  list(mu=mu, Omega=Omega, phi=phi, gamma2=gamma2, t=t, startY=0.5)
}
