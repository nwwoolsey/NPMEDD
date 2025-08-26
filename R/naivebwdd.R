#' Uses SIMEX bandwidth selection in order to fit linear covariates with measurement error to directional responses.
#' Inputs follow the same as npregdd, except hs is a vector of potential bandwidths.
#' @param x covariates
#' @param y responses
#' @param delta values for the regression model to be fitted at
#' @param hs list of potential bandwidths
#' @export
naivebwdd<-function(hs,x,y,delta){
  obj<-function(h,trainx,trainy,testx,testy){ #Objective Function
    r<-naivedd(h,trainx,trainy,testx)
    return(sum(1-cos(testy-r)))
  }
  n<-length(x) 
  remaining<-seq(1,n,1) #List of selected values
  samp1<-sample(remaining,floor(n/5),replace=FALSE) #5 Fold cross validation, selecting test sets from each subset
  samp2<-sample(remaining[-samp1],floor(n/5),replace=FALSE)
  samp3<-sample(remaining[-c(samp1,samp2)],floor(n/5),replace=FALSE)
  samp4<-sample(remaining[-c(samp1,samp2,samp3)],floor(n/5),replace=FALSE)
  samp5<-remaining[-c(samp1,samp2,samp3,samp4)]
  testsamp<-list(samp1,samp2,samp3,samp4,samp5) #Test based on selected subsets
  trainsamp<-list(unlist(testsamp[-1]),unlist(testsamp[-2]),unlist(testsamp[-3]),unlist(testsamp[-4]),unlist(testsamp[-5])) #Train on whats left
  scoreold<-Inf
  for(h in hs){ #Over possible bandwidths
    score<-0
    for(i in 1:5){ #Calculate five fold scores
      score<-obj(h,x[unlist(trainsamp[i])],y[unlist(trainsamp[i])],x[unlist(testsamp[i])],y[unlist(testsamp[i])])+score
      if(is.nan(score)==TRUE){
        score<-Inf
      }
    }
    if(score<scoreold){ #If improvement, save bandwidth
      scoreold<-score
      hout<-h
    }
  }
  r<-naivedd(hout,x,y,delta)
  out<-list("Fit"=r,"h"=hout) #Return fit
  return(out)
}

