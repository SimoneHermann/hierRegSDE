#' Scoring rule of Gneiting and Raftery (...)
#'
#' @description Interval score.
#' @param l lower bound
#' @param u upper bound
#' @param x true value
#' @return interval score
#' @export

scoreRule <- function(l,u,x){
  u-l + 2/0.05*(l-x)*(x<l) + 2/0.05*(x-u)*(x>u)
}
