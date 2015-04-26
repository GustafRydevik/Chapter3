EndemicConstant<-(function(){function(n.infections,incidence,start.time,end.time,...){
  duration=abs(end.time-start.time)
  latest.infection<-rep(NA,n.infections)
  for(i in duration:1){
    latest.infection[sample(1:n.infections,round(n.infections*incidence))]<-i
  }
  return(latest.infection)
}})()



EndemicLinear.old<-(function(){function(n.infections,incidence,start.time,end.time,trend,absolute.ss=TRUE,...){
  ##Absolute.ss indicate the number of positive tests taken?
  duration=abs(end.time-start.time)
  latest.infection<-rep(NA,n.infections)
  ntot<-n.infections
  if(absolute.ss){
    mean.incidence=incidence+(end.time-start.time)*trend/2
    tot.inc=1-(1-mean.incidence)^(end.time-start.time)
    latest.infection<-rep(NA,n.infections/tot.inc)
    ntot<-length(latest.infection)
  }
  incidence.per.time<- (0.5+(0:(duration-1)))*trend+incidence # Trend as measured from Now to Past
  for(i in seq(duration,1)){
    sample.ndx<-sample(1:ntot,round(ntot*incidence.per.time[i]))
    latest.infection[sample.ndx]<-runif(length(sample.ndx),i-1,i)#Ideally a trapezoid..
  }
  if(absolute.ss){
    n.current<-sum(!is.na(latest.infection))
    diff<-n.infections-n.current
    which.na<-which(is.na(latest.infection))
    replace.ndx<-if(diff>0){which.na[seq_len(abs(diff))]}else{
      seq_along(latest.infection)[-which.na][seq_len(abs(diff))]}
    if(diff!=0){latest.infection[replace.ndx]<-sample(latest.infection[-1*sign(diff)*which.na],min(abs(diff),length(which.na)))}
  }
  return(latest.infection)
}})()


EndemicLinear<-(function(){function(n.infections,incidence,start.time,end.time,trend,absolute.ss=TRUE,...){
  ##Absolute.ss indicate the number of positive tests taken?
  duration=abs(end.time-start.time)
  latest.infection<-rep(NA,n.infections)
  ntot<-n.infections
  if(absolute.ss){
    mean.incidence=incidence+(end.time-start.time)*trend/2
    tot.inc=1-(1-mean.incidence)^(end.time-start.time)
    latest.infection<-rep(NA,n.infections/tot.inc)
    ntot<-length(latest.infection)
  }
  
  incidence.per.time<- (0.5+(0:(duration-1)))*trend+incidence # Trend as measured from Now to Past
  for(k in 1:length(latest.infection)){
  for(i in seq(duration,1)){
   infected<-rbinom(1,1,incidence.per.time[i])
   if(infected){
    latest.infection[k]<-runif(1,i-1,i)#Ideally a trapezoid..
  }}}
  if(absolute.ss){
    n.current<-sum(!is.na(latest.infection)) 
    diff<-n.infections-n.current
    which.na<-which(is.na(latest.infection))  ## individuals that have not currently been infected
    
    if(diff>0){  ## if the number of infected individuals is too low
      replace.ndx<-which.na[seq_len(abs(diff))] ## replace some currently uninfected
      latest.infection[replace.ndx]<-sample.int(seq(duration:1),diff,replace=T,prob=incidence.per.time[i])-runif(diff,0,1)
      }
    if(diff<0){ ##otherwise, replace some currently infected
      replace.ndx<-seq_along(latest.infection)[-which.na][seq_len(abs(diff))]
      latest.infection[replace.ndx]<-NA
      }  
  }
  return(latest.infection)
}})()


EndemicLinearWithReinfections<-(function(){function(n.infections,incidence,start.time,end.time,trend,absolute.ss=TRUE,...){
  duration=abs(end.time-start.time)
  latest.infection<-rep(NA,n.infections)
  ntot<-n.infections
  n.inf<-rep(NA,duration)
  incidence.per.time<- (0.5+(0:(duration-1)))*trend+incidence # Trend as measured from Now to Past
  for(i in seq(duration,1)){
    n.inf[i]<-sum(sample(c(0,1),ntot,prob=c(1-incidence.per.time[i],incidence.per.time[i]),replace=T))
  }
  return(n.inf)
}})()



### Demonstration that for Endemiclinear, P(t<T)=exp(a+trend*T/2)
tmp<-hist(EndemicLinear(10000,incidence=0.1,start.time=1,end.time=100,trend=-0.001),
          breaks=0:100,freq=FALSE,ylim=c(0,0.20))$counts
plot(cumsum(tmp)/10000,ylim=c(0,1),type="l")
lines(pexp(1:100,0.1+(1:100*0.001)/2),col="red")


EpidemicExp<-(function(){function(n.infections,start.time=0,end.time,trend=4){
  duration<-end.time-start.time
  scale.factor<-duration
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(trend)-1)*(dist.unif-1/(1-exp(trend))))/trend
  dist.exp<-dist.exp*scale.factor+start.time
  return(dist.exp)
}})()

EpidemicLognorm<-(function(){function(n.infections,end.time,peak.time,sd,...){
  infection.times<-rlnorm(n.infections,meanlog=log(end.time-peak.time),sdlog=log(sd))
  infection.times<-end.time-infection.times
  #infection.times<-exp(infection.times.log)
  return(infection.times)
}})()