end.time<-100
incidence=1/100
pop<-1:1e5
InfTime<-rep(1e5,length(pop))

iteration.data<-
  LabdataGeneratorGeneric(
    testkinetic=kinetic.fun,
    timeFun=EndemicConstant,
    timeFun.args=list(n.infection=1000,
                      start.time=1,
                      end.time=30,
                      incidence=1/100
                      ),
    errorFun=errorFun.lognorm,
    errorFun.args=list(standard.deviation=log(c(1.01,1.01)))#log(measurement.sd))
  )
bugsdata<-list(N=nrow(iteration.data$test.obsvalues),
               Test.data=iteration.data$test.obsvalues,
               time.lookup=seq(0.5,99.5,by=0.5),
               test.lookup=kinetic.fun(c(seq(0.5,99.5,by=0.5))),
               is.naive=is.na(iteration.data$test.obsvalues[,1]),
               ntest=n.tests,censorLimit=30)
pars.inits<-vector(length=n.chains.,mode="list")
for(C in 1:n.chains.){
  pars.inits[[C]]<-list(
    peak.time=rgamma(1,4,scale=Set.endtime/4),
    duration.tmp=abs(rt(1,5)),
    logsd.tmp=rep(runif(1,log(1.02),log(4)),n.tests),
    InfTime=ifelse(is.na(iteration.data$test.obsvalues[,1]),
                   35,runif(nrow(iteration.data$test.obsvalues),0,30))#,
    #alpha=runif(1,0.01,0.3),
    #beta.tmp=runif(1,0,0.01)
  )
}


tmp<-jags.model(file=file.path(script.path,"bugs/hindcast_constant.txt"),
           data=bugsdata,inits=pars.inits,
           n.adapt=100,
           n.chains=5)

tmp.samples<-coda.samples(tmp,variable.names=c("alpha","beta","InfTime","tau"),n.iter=500)
lineardist<-function(time,alpha,beta,chain){
  alpha<-alpha[chain]
  beta<-beta[chain]
  (alpha+beta*time)*exp(-(alpha+beta*time/2)*time)}
