#' @title MC test for AR(1)
#' @name MCtestforAR
#' @description Monte Carlo test for autogressive panel
#' @param data data set
#' @param N number of reputation for the test
#' @param phi the unit roots fot test
#' @param H length of chain of phi of stage 1
#' @param M length of chain of phi of stage 2
#' @return power of value
#' @examples
#' \dontrun{
#' data=ar(0.6,10,rep(0,5))
#' resault=MCtestar(data,0.6,20,20,100)
#' print(resault)
#' }
#' @export
MCtestar=function(data,phi,M,H,N){
  pa=rep(0,N)
  for (i in 1:N) {
    pa[i]=DWH(data,phi,M,H)
  }
  return(mean(pa>0.05))
}
#' @title generate AR(1)
#' @description generate AR(1)
#' @param phi unite root
#' @param t time
#' @param y0 orignal data
#' @return a time series
#' @importFrom  stats rnorm
#' @examples
#' \dontrun{
#' data=ar(0.6,10,rep(0,5))
#' print(data)
#' }
#' @export
ar=function(phi,t,y0){
  m=length(y0)
  A=matrix(0,m,t)
  for (k in 1:m) {
    y=rep(0,t)
    y[1]=(1-phi)*y0[k]+rnorm(1,0,1)
    for (tt in 2:t){
      y[tt]=(1-phi)*y[tt-1]+rnorm(1,0,1)
    }
    A[k,]=y  
  }
  A=cbind(y0,A)
  return(A)
}

pphi=function(data){
  ymean=apply(data[,-c(1,2)],1,mean)
  yh=as.vector(data[,-c(1,2)]-ymean)
  t=dim(data)[2]
  yh.lag=as.vector(data[,-c(1,t)]-ymean)
  phi=1/(yh.lag%*%yh.lag)*(yh.lag%*%yh)
  return(as.numeric( phi))
}

DWH=function(data,phi0,M,H){
  y0=data[,1]
  m=dim(data)[1]
  t=dim(data)[2]-1
  pre.phi=pphi(data)
  cre.phi=rep(0,M+1)
  for (m in 1:(1+M)){
    inst.phi=rep(0,H)
    for (h in 1:H) {
      data0=ar(phi0,t,y0)
      pphi(data0)
      inst.phi[h]=pphi(data0)
    }
    cre.phi[m]=mean(inst.phi)
  }
  Q0=abs(pre.phi-cre.phi[1])
  Qm=abs(cre.phi[2:(M+1)]-cre.phi[1])
  p=mean(Qm>Q0)
  return(p)
}

MCtestar=function(data,phi,M,H,N){
  pa=rep(0,N)
  for (i in 1:N) {
    pa[i]=DWH(data,phi,M,H)
  }
  return(mean(pa>0.05))
}