EndemicDisease<-(function(){function(n.infections,incidence,start.time,end.time,...){
  duration=abs(end.time-start.time)
  latest.infection<-rep(Inf,n.infections)
  for(i in duration:1){
    latest.infection[sample(1:n.infections,round(n.infections*incidence))]<-i
  }
  return(latest.infection)
}})()



EndemicLinear<-(function(){function(n.infections,incidence,start.time,end.time,trend,...){
  duration=abs(end.time-start.time)
  latest.infection<-rep(Inf,n.infections)
  incidence.per.time<- -seq(duration,1)*trend+incidence
  for(i in seq(duration,1)){
    latest.infection[sample(1:n.infections,round(n.infections*incidence.per.time[i]))]<-i
  }
  return(latest.infection)
}})()


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