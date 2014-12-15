
kinetic.fun=LotkaVolterra.Fun(disease=diseaseType3)
C=30
peak=10
obsvalues<-errorFun.lognorm(kinetic.fun(seq(0.01,C,by=0.01)),log(c(1.01,1.01)))

measurement.sd=c(test1.sd=test1.sd,test2.sd=test2.sd)
inftime.bugs<-jags.model(file=file.path(script.path,"bugs/hindcast_only_inftime_test.txt"),
                           data=list(N=nrow(obsvalues),
                                                  time.lookup=c(seq(0.5,99.5,by=0.5)),
                                                  test.lookup=kinetic.fun(c(seq(0.5,99.5,by=0.5))),
                                                  ntest=2,
                                     Test.data=obsvalues,
                                     censorLimit=C,
                                     logsd=log(c(1.01,1.01)),
                                     peak=peak),
                           n.adapt=0,
                           n.chains=1,
                         inits=list(
                           prepeak=sample(c(0,1),nrow(obsvalues),replace=TRUE),
                           Prepeak.time=runif(nrow(obsvalues),0,peak),
                           Postpeak.time=runif(nrow(obsvalues),peak,C)))
                           #InfTime=runif(nrow(obsvalues),0,C)))

inftime.samples<-coda.samples(inftime.bugs,1000,variable.names=c("InfTime"))

inftime.df<-data.frame(true=rep(seq(0.01,C,by=0.01),each=100),est=c(inftime.samples[[1]][901:1000,]))
ggplot(inftime.df,aes(x=factor(round(true)),y=est))+geom_violin()
lv.time.phase
