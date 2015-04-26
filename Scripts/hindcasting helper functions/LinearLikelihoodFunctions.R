#### Reading in  functions specific for this project 


lineardist<-function(time,alpha,beta,chain){
  alpha<-alpha[chain]
  beta<-beta[chain]
  (alpha+beta*time)*exp(-(alpha+beta*time/2)*time)}


dlinear.tsi<-function(InfTime,censorLimit=30,incidence=0.01,trend=0.01/(30*2)){
  area<-(1-exp(-(incidence*censorLimit+trend*censorLimit^2/2)))
  
  pI<-(incidence+trend*InfTime)*exp(-(incidence+trend*(InfTime)/2)*InfTime)/area*(censorLimit>InfTime)
  return(pI)
}


### cumulative probablity of linear.tsi 
rlinear.tsi<-function(n=1000,censorLimit=30,incidence=0.1,trend=0.1/(30*2)){
  if(trend<(-incidence/censorLimit)){stop(simpleError("Too steep linear trend results in \nnegative incidence over the interval"))}
  PI.lookup<-data.frame(I=seq(0,censorLimit,by=0.1),PI=sapply(seq(0,censorLimit,by=0.1),
                                                              function(x){
                                                                integrate(dlinear.tsi,0,x,incidence=incidence,trend=trend)$value}
  ))
  require(data.table)
  runif.dt<-data.table(data.frame(P=runif(n=n)),key="P")
  result.linear.tsi<-data.table(PI.lookup,key="PI")[runif.dt,list(P,I),roll="nearest"]
  return(result.linear.tsi$I)
}

a=0.01
censorLimit=30
b=0.01/30

###Why is the variability much higher in the "direct" generator, than in the simulated one?

par(mfrow=c(2,1))
#hist(round(rlinear.tsi(1000,incidence=a,trend=b)),breaks=0:31)
#hist(EndemicLinear(round(1000/(c*(a+b*c/2))),incidence=a,trend=b,start.time=0,end.time=30),breaks=0:31)

#plot(sapply(seq(0,30,by=0.1),function(x)integrate(dlinear.tsi,x,30)$value))
#test.data1<-rlinear.tsi(1000,incidence=a,trend=b,censorLimit=censorLimit)
#test.data2<-EndemicLinear(round(1000/(censorLimit*(a+b*censorLimit/2))),incidence=a,trend=-b,start.time=0,end.time=censorLimit)


#A function for estimating the RMSEP between a the linear trend curve decided by pars, and the histogram of data
tsi.linear.rmse<-function(pars,censorLimit.=censorLimit,data=obs.data,transform.pars=TRUE){
  incidence.<-pars[1]
  trend.<-pars[2]
  if(transform.pars){
    incidence.<-exp(incidence.)
    trend.<-exp(trend.)-incidence./censorLimit.
  }
  data.prop<-hist(data,breaks=0:censorLimit.,plot=F)
  predprop<-dlinear.tsi(data.prop$mids,incidence=incidence.,trend=trend.,censorLimit=censorLimit.)
  predprop.se<-(predprop-data.prop$density)^2
  return(sqrt(mean(predprop.se)))
}

#A function for evaluation the data likelihood of the linear trend curve decided by pars, given data.
#Accepts parameters>0, ie logtransformed, so as to work well in conjunction with nlm
##Popsize gives the number of total individuals tested to generate the observed number of cases. 


tsi.linear.ll<-function(pars,censorLimit.,data=obs.data,transform.pars=TRUE,popsize,return.components=FALSE){
  incidence.<-pars[1]
  trend.<-pars[2]
  if(transform.pars){
    incidence.<-exp(incidence.)
    trend.<-exp(trend.)-incidence./censorLimit.
  }
  InfTime<-data[!is.na(data)]
  mean.incidence<-incidence.+trend.*censorLimit./2  ##calculating the average incidence over the period
  outside.range<-any(incidence.+trend.*InfTime<0)  ## checks if the estimated incidence is negative
  
  ##Calculates the evaluation of the prior for the trend line
  trend.prior.ll<-dunif(trend.,-2*mean.incidence/censorLimit.,2*mean.incidence/censorLimit.,log=T)
 
  ##Calculates the value of the beta prior for the mean incidence over the period... (is this the corrected version?) 
  mean.inc.ll<-dbeta(1-(1-mean.incidence)^censorLimit.,length(InfTime)+1,(popsize-length(InfTime))+1,log=T)
  
  ##calulating the normalizing area factor. 
  area<-log((1-exp(-(incidence.+trend.*censorLimit./2)*censorLimit.)))
  
  
  ##calculating the contribution of the cdf of the time-since-infection distribution
  linear.tsi.ll<-sum(log((incidence.+trend.*InfTime)*exp(-(incidence.+trend.*(InfTime)/2)*InfTime))-area)
  
  loglik<-trend.prior.ll+mean.inc.ll+
    linear.tsi.ll+
    log((1-outside.range))
  ret.value<-loglik
  if(return.components){
    ret.value<-list(area=area,
                    trend.prior.ll=trend.prior.ll,
                    mean.inc.ll=mean.inc.ll,
                    tsi.ll=linear.tsi.ll-area)
    
  }
  return(ret.value)
}


tsi.beta.lambda.ll<-function(pars,censorLimit.,data=obs.data,transform.pars=TRUE,popsize,return.components=FALSE){
  mean.incidence<-pars[1]
  trend.<-pars[2]
  if(transform.pars){
    mean.incidence<-exp(mean.incidence)
    trend.<-exp(trend.)-mean.incidence/censorLimit.*2
  }
  InfTime<-data[!is.na(data)]
  #mean.incidence<-incidence.+trend.*censorLimit./2
  incidence.<-mean.incidence-trend.*censorLimit./2
  outside.range<-any(incidence.+trend.*InfTime<0)
  trend.prior.ll<-dunif(trend.,-2*mean.incidence/censorLimit.,2*mean.incidence/censorLimit.,log=T)
  mean.inc.ll<-dbeta(1-(1-mean.incidence)^censorLimit.,length(InfTime)+1,(popsize-length(InfTime))+1,log=T)
  area<-log((1-exp(-(incidence.+trend.*censorLimit./2)*censorLimit.)))
  linear.tsi.ll<-sum(log((incidence.+trend.*InfTime)*exp(-(incidence.+trend.*(InfTime)/2)*InfTime))-area)
  
    loglik<-trend.prior.ll+mean.inc.ll+
   +log((1-outside.range))
  ret.value<-loglik
  if(return.components){
    ret.value<-list(area=area,
                    trend.prior.ll=trend.prior.ll,
                    mean.inc.ll=mean.inc.ll,
                    tsi.ll=linear.tsi.ll-area)
  }
  return(ret.value)
}


#A function for attempting to use nlm to fit the linear tsi curve. 
tsi.linear.est<-function(incidence=0.01,
                         trend=0.01/30,
                         censorLimit=30,
                         simtype="population",
                         n.obs=1000,engine="nlm",transform.pars=TRUE,type="ll"){
  
  
  approx.popsize<-round(n.obs/(censorLimit*(incidence+trend*censorLimit/2)))
  obs.data<-if(simtype=="population"){EndemicLinear(approx.popsize,incidence=incidence,trend=trend,start.time=0,end.time=censorLimit)
  }else{rlinear.tsi(n.obs,incidence=incidence,trend=trend,censorLimit=censorLimit)}
  
  obs.data<-obs.data[!is.na(obs.data)]
  prior.inc<-exp(censorLimit*log(1-length(obs.data)/(approx.popsize)))
  pars.est<-nlm(tsi.linear.ll,p=c(log(prior.inc),
                                  log((prior.inc)/(censorLimit))),
                transform.pars=transform.pars,data=obs.data,popsize=approx.popsize,censorLimit.=censorLimit,
                steptol=1e-8)
  if(transform.pars){
    pars.est$estimate[1]<-exp(pars.est$estimate[1])
    pars.est$estimate[2]<-exp(pars.est$estimate[2])-pars.est$estimate[1]/censorLimit
  }
  return(list(incidence.est=pars.est$estimate[1],
              trend.est=pars.est$estimate[2],
              incidence.bias=incidence-pars.est$estimate[1],
              trend.bias=trend-pars.est$estimate[2])
  )
}

#lambda.hat<-a+b.abs*c/2
#tsi.beta.lambda.ll(pars=c(log(lambda.hat),log(b.abs+lambda.hat*c/2)),censorLimit. = c,data = plot.data,transform.pars = TRUE,popsize = 7000)

