#' Default values
#'
#' @description Default values for \eqn{\mu} for the specific model out of {Gompertz, Richards, logistic, Weibull}
#' @param model one out of {SDE, ODE}
#' @param mod one out of {Gompertz, Richards, logistic, Weibull}
#' @param para prior values for the data set of Virkler et al. (1979) or parameters for the simultions
#' @return vector of parameters
#' @export

getPar <- function(model="ODE", mod="Gompertz", para=c("prior","truePar")){
  if(para=="prior"){

    if(model=="ODE"){

      if(mod=="Gompertz"){
        m <- c(30,5,7)
      }else{
        if(mod=="logistic"){
          m <- c(30,9,7)
        }else{
          if(mod=="Richards"){
            m <- c(25,0.1,5,40)
          }else{
            if(mod=="Weibull"){
              m <- c(28,20,8,5)
            }else{
              if(mod=="Paris"){
                m <- c(25,2,1.8)
              }else{
                if(mod=="Paris2"){
                  m <- c(50,-10)
                }
              }
            }
          }
        }
      }
    }else{ # SDE
      if(mod=="Gompertz"){
        m <- c(5,7)
      }else{
        if(mod=="logistic"){
          m <- c(30,7)
        }else{
          if(mod=="Richards"){
            m <- c(0.1,5,40)
          }else{
            if(mod=="Weibull"){
              m <- c(28,8,5)
            }else{
              if(mod=="Paris"){
                m <- c(25,1.8)
              }else{
                if(mod=="Paris2"){
                  m <- c(-10,1)  # codes are only for vectors... for diagonal Omega, it does not change anything
                }
              }
            }
          }
        }
      }
    }

  }else{  # for simulation

    if(model=="ODE"){

      if(mod=="Gompertz"){
        m <- c(10,3,5)
      }else{
        if(mod=="logistic"){
          m <- c(10,9,7)
        }else{
          if(mod=="Richards"){
            m <- c(10,2,8,3)
          }else{
            if(mod=="Weibull"){
              m <- c(10,9,5,2)
            }else{
              if(mod=="Paris"){
                m <- c(10,1,1.2)
              }else{
                if(mod=="Paris2"){
                  m <- c(20,-3)
                }
              }
            }
          }
        }
      }
    }else{  # SDE

      if(mod=="Gompertz"){
        m <- c(3,5)
      }else{
        if(mod=="logistic"){
          m <- c(10,7)
        }else{
          if(mod=="Richards"){
            m <- c(2,8,3)
          }else{
            if(mod=="Weibull"){
              m <- c(10,5,2)
            }else{
              if(mod=="Paris"){
                m <- c(10,1.2)
              }else{
                if(mod=="Paris2"){
                  m <- c(20,-3)  # codes are only for vectors... for diagonal Omega, it does not change anything
                }
              }
            }
          }
        }
      }
    }

  }
  m
}
