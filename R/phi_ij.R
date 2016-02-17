#' Helping function
#'
#' @description transfers the list to one vector
#' @param phi list of samples for mixed effect
#' @param i component of vector
#' @param j which individual
#'
#' @return samples for the i-th component of \eqn{\phi_j}
#' @export



phi_ij <- function(phi,i,j){
  sapply(phi, function(mat){mat[i,j]})
}
