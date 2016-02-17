#' Variance function
#'
#' @description Gives out the variance function (SDE)
#' @param modelVar one out of {, AR}
#' @return function
#' @export
getFunVar <- function(modelVar=""){
  sVar <- function(t,x=1) 1
  if(modelVar=="AR") sVar <- function(t, x) x
  sVar
}
