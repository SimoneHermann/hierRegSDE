#' Posterior for mu
#'
#' @description Posterior for parameters \eqn{\mu}
#' @param phi
#' @param mu
#' @param v
#' @param Omega
#' @return one sample of posterior



postmu <- function(phi, m, v, Omega){  # phi matrix
  n <- nrow(phi)
  V_ <- diag(1/v)
  Vpost <- solve(V_ + n*solve(Omega))
  mpost <- Vpost%*%( V_%*%m + apply((solve(Omega)%*%t(phi)),1,sum) )

  rmvnorm(1,mpost,Vpost)
}
