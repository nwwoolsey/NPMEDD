#Example for Nonparametric Regression for A Circular Response with Error-in-Covariate
#This will take a long time to compile (on the order of 3-20 hours). It can be sped up by reducing the width of potential bandwidths, and by reducing the
#replicates for the complex error method, but that will also increase how wiggly the fit is.
set.seed(1)
library(lpme)
l<-read.csv("Data/tnrcc_wind_direction_resultant.txt")
convert<-function(x){
  o<-as.POSIXct(l[,2],format="%H:%M")
  out<-as.integer(format(o,format="%H"))
  return(out)
}
time<-convert(l[,2])
l<-l[,-(1:2)]
direction<-matrix(l,ncol=1,nrow=1756*123)
direction<-as.integer(l[,80])*pi/180
timenew<-time[-which(is.na((direction)))]
direction<-direction[-which(is.na((direction)))]
sigmau<-2.306
timeerr<-timenew+rnorm(length(timenew),0,sigmau)
g<-seq(3,22,.5)
tru<-naivebwdd(seq(.5,5,.1),timenew,direction-pi,g)
naive<-naivebwdd(seq(1,5,.1),timeerr,direction-pi,g)
dec<-npddbw(timeerr,direction-pi,seq(.55,.65,.01),"deconv",g,"normal",sigmau,500)
compbw<-npddbw(timeerr,direction-pi,seq(.92,.97,.01),"comp",g,"normal",sigmau,500)$h#we use a reduced number of replicates to find the optimal bandwidth to save computation time
comp<-npregdd(compbw,"comp",timeerr,direction-pi,g,"normal",sigmau,10000)
os<-npddbw(timeerr,direction-pi,seq(.6,.75,.01),"onestep",g,"normal",sigmau,500)
plot(timeerr,direction-pi,xlim=c(0,24),ylim=c(-1.7,.8),xlab="Time of Day",ylab="Wind Direction In Radians",cex=.35)
lines(g,tru$Fit)
lines(g,dec$Fit,lty=2,col="blue")
points(g,dec$Fit,pch=1,col="blue",cex=.75)
lines(g,comp,lty=3,col="orange")
points(g,comp,pch=2,cex=.75,col="orange")
lines(g,os$Fit,lty=4,col="purple")
points(g,os$Fit,pch=5,col="purple",cex=.75)
lines(g,naive$Fit,lty=5,col="red")
points(g,naive$Fit,pch=0,col="red",cex=.75)
legend("topright",col=c("black","blue","orange","purple","red"),legend=c("No Error","Deconvoluting","Complex Error Corrrection","One Step","Naive"),lty=c(1,2,3,4,5),cex=1)

