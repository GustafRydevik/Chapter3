EndemicDisease<-(function(){function(n.infections,incidence,start.time,end.time,...){
  duration=abs(end.time-start.time)
  latest.infection<-rep(NA,n.infections)
  for(i in duration:1){
    latest.infection[sample(1:n.infections,round(n.infections*incidence))]<-i
  }
  return(latest.infection)
}})()



EndemicLinear<-(function(){function(n.infections,incidence,start.time,end.time,trend,...){
  duration=abs(end.time-start.time)
  latest.infection<-rep(NA,n.infections)
  incidence.per.time<- -(1:duration)*trend+incidence
  for(i in seq(duration,1)){
    latest.infection[sample(1:n.infections,round(n.infections*incidence.per.time[i]))]<-i
  }
  return(latest.infection)
}})()


### Demonstration that for Endemiclinear, P(t<T)=exp(a+trend*T/2)
tmp<-hist(EndemicLinear(10000,incidence=0.1,start.time=1,end.time=100,trend=-0.001),
          breaks=0:100,freq=FALSE,ylim=c(0,0.20))$counts
plot(cumsum(tmp)/10000,ylim=c(0,1),type="l")
lines(pexp(1:100,0.1+(1:100*0.001)/2),col="red")


EpidemicExp<-(function(){function(n.infections,start.time=0,end.time,slope=4){
  duration<-end.time-start.time
  scale.factor<-duration
  dist.unif<-runif(n.infections,0,1)
  dist.exp<-log((exp(slope)-1)*(dist.unif-1/(1-exp(slope))))/slope
  dist.exp<-dist.exp*scale.factor+start.time
  return(dist.exp)
}})()



EpidemicLognorm<-(function(){function(n.infections,end.time,peak.time,sd){
  infection.times<-rlnorm(n.infections,meanlog=log(end.time-peak.time),sdlog=log(sd))
  infection.times<-end.time-infection.times
  #infection.times<-exp(infection.times.log)
  return(infection.times)
}})()