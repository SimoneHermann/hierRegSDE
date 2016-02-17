#' Example function for the nonlinear regression model
#'
#' @description Example function - how to use the functions in the package
#' @param mod model out of {Gompertz, Richards, logistic, Weibull}
#' @param mixed default=TRUE; if false: \eqn{\phi_i=\mu}
#' @return nice pictures
#' @export

testReg <- function(mod="Gompertz", mixed=TRUE){
  fODE <- getFun("ODE", mod)
  mu <- getPar("ODE", mod, "truePar")

  n <- 5
  parameters <- defaultPar(mu, n, mixed)
  y <- drawData("ODE", fODE, parameters)
  t <- parameters$t
  lt <- length(t)

  prior <- getPrior(mu, parameters$gamma2, mixed)
  start <- getStart(mu, n, mixed)
  cut <- ceiling(lt/2)
  thinned <- seq(101,1100,by=2)

  if(mixed){
    Est_ode <- estReg(t, y, prior=prior, start=start, mod=mod, ipred=1, cut=cut, len=1100)
    samples <- list(phi=t(sapply(Est_ode$phi[thinned], function(mat) mat[1,])), gamma2=Est_ode$gamma2[thinned])
  }else{
    Est_ode <- estReg_single(t, y[1,], prior=prior, start=start, mod=mod, len=1100)
    samples <- list(phi=Est_ode$phi, gamma2=Est_ode$gamma2)
  }
  ypred <- predReg(t[(cut+1):lt], samples, mod=mod, cand=y[1,cut]+seq(-2,5,by=0.01), len=100)
  qu1 <- apply(ypred,2,quantile,0.025)
  qu2 <- apply(ypred,2,quantile,0.975)
  mit <- apply(ypred,2,mean)
  isIn <- y[1,(cut+1):lt]<=qu2 & y[1,(cut+1):lt]>=qu1
  par(mfrow=c(2,2))
  if(mixed){
    plot(Est_ode$mu[,1],ylim=range(mu)+c(-2,2),main="Markov Chains",ylab=expression(mu))
    for(i in 2:length(mu)){
      points(Est_ode$mu[,i])
    }
    abline(h=mu,col=2)
    acf(phi_ij(Est_ode$phi,2,1),main=expression("Acf of simulations for "~phi[1]~" in the second series"),cex=0.5)
  }else{
    plot(Est_ode$phi[,1],ylim=range(mu)+c(-2,2),main="Markov Chains",ylab=expression(phi))
    for(i in 2:length(mu)){
      points(Est_ode$phi[,i])
    }
    abline(h=mu,col=2)
    acf(Est_ode$phi[,1],main=expression("Acf of simulations for "~phi[1]),cex=0.5)
  }
  hist(Est_ode$gamma2, main=expression("Posterior of "~sigma^2),xlab=expression(sigma^2),freq=F,breaks=seq(min(Est_ode$gamma2),max(Est_ode$gamma2),length=50))
  abline(v=parameters$gamma2, col=2)
  plot(t,y[1,],pch=20,ylab="y",ylim=c(min(y),max(qu2)))
  for(i in 2:n){
    points(t,y[i,],col="grey")
  }
  points(t,y[1,],pch=20)
  lines(t[(cut+1):lt],qu1,col=2)
  lines(t[(cut+1):lt],qu2,col=2)
  lines(t[(cut+1):lt],mit,col=4)
  legend("bottomright",c("predicted series","Quantiles","Mean"),col=c(1,2,4),pch=c(20,-1,-1),lty=c(-1,1,1),cex=0.7,box.lty=0)
  print(paste("amount of covering intervals:", mean(isIn)))
}
