#' Yields the regression or drift function for the specific models
#'
#' @description Gives out the regression function (ODE) or the drift function (SDE) for the specific model out of {Gompertz, Richards, logistic, Weibull}
#' @param model one out of {SDE, ODE}
#' @param mod one out of {Gompertz, Richards, logistic, Weibull}
#' @return function
#' @export

getFun <- function(model, mod){
  if(mod=="Gompertz"){
    fODE <- function(phi,t){
      phi[1]*exp(-phi[2]*exp(-phi[3]*t))
    }
    bSDE <- function(phi,t,x){
      phi[1]*phi[2]*exp(-phi[2]*t)*x
    }
  }else{
    if(mod=="logistic"){
      fODE <- function(phi,t){
        phi[1]/(1+phi[2]*exp(-phi[3]*t))
      }
      bSDE <- function(phi,t,x){
        phi[2]*x*(1-1/phi[1]*x)
      }
    }else{
      if(mod=="Richards"){
        fODE <- function(phi,t){
          phi[1]/(1+phi[2]*exp(-phi[3]*t))^phi[4]
        }
        bSDE <- function(phi,t,x){
          phi[1]*phi[2]*phi[3]*exp(-phi[2]*t)/(1+phi[1]*exp(-phi[2]*t))*x
        }
      }else{
        if(mod=="Weibull"){
          fODE <- function(phi,t){
            phi[1]-phi[2]/exp(phi[3]*t^phi[4])
          }
          bSDE <- function(phi,t,x){
            phi[2]*phi[3]*t^(phi[3]-1)*(phi[1]-x)
          }
        }else{
          if(mod=="Paris"){
            fODE <- function(phi,t){
              phi[1]-phi[2]*t^(-phi[3])
            }
            bSDE <- function(phi,t,x){
              phi[2]/t*(phi[1]-x)
            }
          }else{
            if(mod=="Paris2"){
              fODE <- function(phi,t){
                phi[1]*log(t)+phi[2]
              }
              bSDE <- function(phi,t,x){
                (x-phi[1])/(t*log(t))
              }
            }
          }
        }
      }
    }
  }
  if(model=="ODE"){
    return(fODE)
  }else{
    return(bSDE)
  }

}
