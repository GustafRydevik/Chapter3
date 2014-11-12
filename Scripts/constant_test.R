#proj.path<-"/Users/gustafrydevik/Dropbox/PhD folder/Chapter3"

end.time<-30
incidence<-1/10
Pop<-1:1000
InfTime<-rep(NA,length(Pop))

kinetic.fun=LotkaVolterra.Fun(diseaseType1)
iteration.data<-
  LabdataGeneratorGeneric(
    testkinetic=kinetic.fun,
    timeFun=EndemicConstant,
    timeFun.args=list(n.infection=1000,
                      start.time=1,
                      end.time=30,
                      incidence=1/10
    ),
    errorFun=errorFun.lognorm,
    errorFun.args=list(standard.deviation=log(c(1.2,1.1)))#log(measurement.sd))
  )


library(rjags)
tmp<-jags.model(file.path(script.path,"bugs/hindcast_constant.txt"),
           data=list(N=nrow(iteration.data$test.obsvalues),
                     Test.data=iteration.data$test.obsvalue,
             is.naive=as.numeric(is.na(iteration.data$infection.times)),
             censorLimit=end.time,ntest=2,
             time.lookup=seq(0.5,99.5,by=0.5),
             test.lookup=kinetic.fun(c(seq(0.5,99.5,by=0.5)))),
           inits=list(InfTime=ifelse(is.na(iteration.data$infection.times),end.time+1,NA),lambda=1/10),
                n.chains=5,n.adapt=100)


tmp2<-coda.samples(tmp,c("InfTime","incidence","sd"),n.iter=1000,n.adapt=500)


#diagnostics
plot(tmp2[,c("incidence","sd[1]","sd[2]")])
hist(unlist(tmp2[[1]][,grep("InfTime",colnames(tmp2[[1]]))]),freq=FALSE)
lines(dexp(seq(0,end.time,by=0.5),mean(tmp2[[1]][,"incidence"])))
plot(iteration.data$infection.times,colMeans(tmp2[[1]][,grep("InfTime",colnames(tmp2[[1]]))]))
