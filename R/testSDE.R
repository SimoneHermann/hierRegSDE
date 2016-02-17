#' Example function for the SDE model
#'
#' @description Example function - how to use the functions in the package
#' @param mod model out of {Gompertz, Richards, logistic, Weibull}
#' @param mixed default=TRUE; if false: \eqn{\phi_i=\mu}
#' @return nice pictures
#' @export

testSDE <- function(mod="Gompertz", mixed=TRUE){
  bSDE <- getFun("SDE", mod)
  mu <- getPar("SDE", mod, "truePar")

  n <- 5
  parameters <- defaultPar(mu, n, mixed)
  Y <- drawData("SDE", bSDE, parameters)
  t <- parameters$t
  lt <- length(t)

  prior <- getPrior(mu, parameters$gamma2, mixed)
  start <- getStart(mu, n, mixed)
  cut <- ceiling(lt/2)
  thinned <- seq(101,1100,by=2)
  if(mixed){
    Est_sde <- estSDE(t, Y, mod=mod, ipred=1, cut=cut, prior=prior, start=start, len=1100)
    samples <- list(phi=t(sapply(Est_sde$phi[thinned], function(mat) mat[1,])), gamma2=Est_sde$gamma2[thinned])
  }else{
    Est_sde <- estSDE_single(t, Y[1,], mod=mod, prior=prior, start=start, len=1100)
    samples <- list(phi=Est_sde$phi[thinned,], gamma2=Est_sde$gamma2[thinned])
  }
  pred <- predSDE(t[cut:lt], samples, last=Y[1,cut], mod=mod)
  qu1 <- apply(pred,1,quantile,0.025)
  qu2 <- apply(pred,1,quantile,0.975)
  mit <- apply(pred,1,mean)
  isIn <- Y[1,(cut+1):lt]<=qu2 & Y[1,(cut+1):lt]>=qu1
  par(mfrow=c(2,2))
  if(mixed){
    plot(Est_sde$mu[,1],ylim=range(mu)+c(-1,1),main="Markov Chains",ylab=expression(mu))
    for(i in 2:length(mu)){
      points(Est_sde$mu[,i])
    }
    abline(h=mu,col=2)
    acf(phi_ij(Est_sde$phi,2,1),main=expression("Acf of simulations for "~phi[1]~" in the second series"),cex=0.5)
  }else{
    plot(Est_sde$phi[,1],ylim=range(mu)+c(-1,1),main="Markov Chains",ylab=expression(phi))
    for(i in 2:length(mu)){
      points(Est_sde$phi[,i])
    }
    abline(h=mu,col=2)
    acf(Est_sde$phi[,1],main=expression("Acf of simulations for "~phi[1]),cex=0.5)
  }
  hist(Est_sde$gamma2, main=expression("Posterior of "~gamma^2),xlab=expression(gamma^2),freq=F,breaks=seq(min(Est_sde$gamma2),max(Est_sde$gamma2),length=50))
  abline(v=parameters$gamma2, col=2)
  plot(t,Y[1,],pch=20,ylab="X",ylim=c(min(Y),max(max(Y),max(qu2))))
  for(i in 2:n){
    points(t,Y[i,],col="grey")
  }
  points(t,Y[1,],pch=20)
  lines(t[(cut+1):lt],qu1,col=2)
  lines(t[(cut+1):lt],qu2,col=2)
  lines(t[(cut+1):lt],mit,col=4)
  legend("bottomright",c("predicted series","Quantiles","Mean"),col=c(1,2,4),pch=c(20,-1,-1),lty=c(-1,1,1),cex=0.7,box.lty=0)
  print(paste("amount of covering intervals:", mean(isIn)))
}
