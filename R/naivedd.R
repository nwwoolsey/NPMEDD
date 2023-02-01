#' Fit linear covariates to directional responses nonparametrically.
#' h is a chosen bandwidth, x is a linear covariate, y is a directional response, and delta is a vector of points to fit the model over.
#' @export
naivedd<-function(h,x,y,delta){#non corrected weight function
  r<-c()
  for(u in delta){
    K<-function(l){
      return(1/sqrt(2*pi)*exp(-l^2/2))#std normal
    }
    g1<-mean(sum(sin(y)*K((u-x)/h)/h)*(1/h*sum(K((u-x)/h)*(u-x)^2)-(u-x)*sum(1/h*K((u-x)/h)*(u-x))))
    g2<-mean(sum(cos(y)*K((u-x)/h)/h)*(1/h*sum(K((u-x)/h)*(u-x)^2)-(u-x)*sum(1/h*K((u-x)/h)*(u-x))))
    m<-atan(g1/g2)
    r<-append(r,m)
  }
  return(r)
}
