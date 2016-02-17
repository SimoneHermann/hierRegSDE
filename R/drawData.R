#' Function for simulating the data
#'
#' @description Simulation of variables depend on the model.
#' @param model one out of {SDE, ODE}
#' @param fun regression or drift function
#' @param para list of parameters, for example defaultPar(getPar(model, para=truePar))
#' @param sVar variance function
#' @param mw mesh width (in the case of SDE)
#' @return data series in para$t
#' @export

drawData <- function(model, fun, para, sVar, mw=10){
  if(!(model%in%c("ODE", "SDE"))){
    print("choose model ODE or SDE")
    return(NULL)
  }
  lt <- length(para$t)
  n <- nrow(para$phi)
  if(missing(sVar)) sVar <- function(t,x=0) 1
  if(model=="ODE"){
    y <- matrix(0,n,lt)
    for(i in 1:n){
      y[i,] <- fun(para$phi[i,],para$t) + rnorm(lt,0,sqrt(para$gamma2)*sVar(para$t))
    }
  }
  if(model=="SDE"){
    lt2 <- (lt-1)*mw+1
    t <- seq(min(para$t),max(para$t),length=lt2)
    dt <- t[2]-t[1]
    Y <- matrix(0,n,lt2)
    vec <- numeric(lt2)

    for(i in 1:n){
      err <- rnorm(lt2-1,0,sqrt(dt))
      vec[1] <- para$startY
      for(j in 2:lt2){
        vec[j] <- vec[j-1] + fun(para$phi[i,],t[j-1],vec[j-1])*dt + sqrt(para$gamma2)*sVar(t[j-1],vec[j-1])*err[j-1]
      }
      Y[i,] <- vec
    }
    ind <- seq(1,lt2,by=mw)
    y <- Y[,ind]
  }
  y
}
