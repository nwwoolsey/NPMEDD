#' Uses SIMEX bandwidth selection in order to fit linear covariates with measurement error to directional responses.
#' Inputs follow the same as npregdd, except hs is a vector of potential bandwidths.
#' @param x covariates
#' @param y responses
#' @param hs list of potential bandwidths
#' @param method "comp" for complex errors "onestep" for one step estimator or "deconv" for deconvoluting estiamtor
#' @param delta values for the regression model to be fitted at
#' @param error error distribution, "normal" or "laplace"
#' @param rep number of repetitions for monte carlo convergence of complex method
#' @param sigmau standard deviation of the measurement error
#' @export
npddbw<-function(x,y,hs,method,delta,error,sigmau,rep){
  if(missing(rep)){
    rep<-200
  }
  if(missing(error)){
    error<-"normal"
  }
  if(sigmau<=0){
    return("Measurement Error Must be Positive")
  }
  n<-length(x)
  CV1<-c()
  CV2<-c()
  err1<-stats::rnorm(n,mean=0,sd=sigmau)
  err2<-err1+stats::rnorm(n,mean=0,sd=sigmau)
  wstar<-x+err1
  wstar2<-x+err2
  for(h in hs){
    cv1<-0
    cv2<-0
    remaining<-seq(1,n,1)
    samp1<-sample(remaining,floor(n/5),replace=FALSE)
    samp2<-sample(remaining[-samp1],floor(n/5),replace=FALSE)
    samp3<-sample(remaining[-c(samp1,samp2)],floor(n/5),replace=FALSE)
    samp4<-sample(remaining[-c(samp1,samp2,samp3)],floor(n/5),replace=FALSE)
    samp5<-remaining[-c(samp1,samp2,samp3,samp4)]
    testsamp<-list(samp1,samp2,samp3,samp4,samp5)
    trainsamp<-list(unlist(testsamp[-1]),unlist(testsamp[-2]),unlist(testsamp[-3]),unlist(testsamp[-4]),unlist(testsamp[-5]))
    for(i in 1:5){
      fit1<-npregdd(h,method,wstar[unlist(trainsamp[i])],y[unlist(trainsamp[i])],x[unlist(testsamp[i])],error,sigmau,rep)
      cv1<-sum(1-cos(y[unlist(testsamp[i])]-fit1))+cv1
      fit2<-npregdd(h,method,wstar2[unlist(trainsamp[i])],y[unlist(trainsamp[i])],x[unlist(testsamp[i])],error,sigmau,rep)
      cv2<-sum(1-cos(y[unlist(testsamp[i])]-fit2))+cv2
    }
    CV1<-append(CV1,cv1)
    CV2<-append(CV2,cv2)
  }
  h1<-hs[which(CV1==min(stats::na.omit(CV1)))]
  h2<-hs[which(CV2==min(stats::na.omit(CV2)))]
  out<-h1^2/h2
  hopt<-hs[which(abs(hs-out)==min(abs(hs-out)))]
  r<-npregdd(h,method,x,y,delta,error,sigmau,rep)
  return(list("h"=hopt,"Fit"=r))
}



