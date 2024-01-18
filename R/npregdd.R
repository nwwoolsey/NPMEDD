#' Fits linear covariates with measurement error to directional responses, has three different methods.
#' h is bandwidth, x is linear covariates, y is directional responses, delta is a vector to fit the model over, sigmau is the standard deviation of the measurement error.
#' This function has three potential methods, "comp" for complex error correction, "onestep" for onestep error correction, and "deconv" for error correction using deconvoluting kernels.
#' rep is the number of random normal imaginary errors to generate per replication for complex error correction, increasing this will smooth the fit, but will increase computation time.
#' @param x covariates
#' @param y responses
#' @param h bandwidth
#' @param method "comp" for complex errors "onestep" for one step estimator or "deconv" for deconvoluting estiamtor
#' @param delta values for the regression model to be fitted at
#' @param error error distribution, "normal" or "laplace"
#' @param rep number of repetitions for monte carlo convergence of complex method
#' @param sigmau standard deviation of the measurement error
#' @export
npregdd<-function(h,method,x,y,delta,error,sigmau,rep){
  if(missing(rep)){
    rep<-200
  }
  if(missing(error)){
    error<-"normal"
  }
  if(sigmau<=0){
    return("Measurement Error Must be Positive")
  }
  if(h<=0){
    return("Bandwidth Must be Positive")
  }
  if(method=="deconv"){
    g1<-sin(y)
    g2<-cos(y)
    m1<-lpme::meanreg(g1,x,h,xgrid=delta,method="DFC",sig=sigmau,error=error)$yhat
    m2<-lpme::meanreg(g2,x,h,xgrid=delta,method="DFC",sig=sigmau,error=error)$yhat
    m<-atan(m1/m2)
    return(m)
  }
  if(method=="comp"){
    r<-c()
    onea<-rep(1,length(x))
    oneb<-rep(1,length(delta))
    diff<-(delta%*%t(onea)-oneb%*%t(x))
    lol<-function(x){
      return(1/sqrt(2*pi)*exp(-(x^2)/(2^2)))
    }
    K<-function(s){
      kern<-(lol(s/h)*(sum(lol(s/h)*s^2)-s*sum(lol(s/h)*s)))/(sum(lol(s/h)*s^2)*sum(lol(s/h))-sum(lol(s/h)*s)^2)
      return(kern)
    }
    errs<-function(s,h){
      t<-matrix(s,ncol=rep,nrow=length(y))+1i*sigmau*matrix(stats::rnorm(length(y)*rep,mean=0,sd=1),ncol=rep,nrow=length(y))
      out<-Re(apply(apply(t,2,K),1,mean))
      return(out)
    }
    kernels<-apply(diff,1,errs)
    g1<-sin(y)%*%kernels
    g2<-cos(y)%*%kernels
    m<-atan(g1/g2)
    r<-append(r,m)
    return(r)
  }
  if(method=="onestep"){
    g1<-sin(y)
    g2<-cos(y)
    m1<-lpme::meanreg(g1,x,h,xgrid=delta,method="HZ",sig=sigmau,error=error)$yhat
    m2<-lpme::meanreg(g2,x,h,xgrid=delta,method="HZ",sig=sigmau,error=error)$yhat
    m<-atan(m1/m2)
    return(m)
  }
}
