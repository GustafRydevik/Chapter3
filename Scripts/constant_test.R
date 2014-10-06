proj.path<-"/Users/gustafrydevik/Dropbox/PhD folder/Chapter3"

end.time<-30
incidence<-1/10
Pop<-1:1000
InfTime<-rep(NA,length(Pop))
for(i in end.time:1){InfTime[sample(Pop,rbinom(1,length(Pop),incidence),replace=F)]<-i}
library(rjags)
tmp<-jags.model(file.path(proj.path,"Scripts/bugs/hindcast_constant_test.txt"),
           data=list(InfTime=InfTime,N=length(InfTime),is.censored=as.numeric(is.na(InfTime)),censorLimit=end.time),
           inits=list(InfTime=ifelse(is.na(InfTime),end.time+1,NA),lambda=1/10))
tmp2<-jags.samples(tmp,c("lambda"),n.iter=1000)
plot(tmp2[[1]][,,],type="l")
