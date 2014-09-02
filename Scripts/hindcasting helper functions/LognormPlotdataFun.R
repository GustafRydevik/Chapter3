lognorm.plotdata.fun<-function(mcmc.posterior,interval.eval=1,range.eval=1000,df=TRUE,use.z=F,...){
  posterior.allchains<-do.call("rbind",mcmc.posterior)
  if(!use.z){
    lognorm.density.matrix<-apply(posterior.allchains,1,function(x){
      dlnorm(seq(0,range.eval,by=interval.eval),meanlog=log(x["peak.time"]),sdlog=log(x["duration"]))})}
  if(use.z){
    X<-seq(0,range.eval,by=interval.eval)    
    logMu<-rep(log(posterior.allchains[,"peak.time"]),each=length(X)) ##Keep each value constant for the range of X
    logSD<-rep(log(posterior.allchains[,"duration"]/2),each=length(X))
    
    X.rep<-rep(X,nrow(posterior.allchains))
    
    Z<-(X.rep/exp(logMu))^(1/logSD)
    scaling.log<-(log(X.rep)-logMu)/logSD-log(X.rep)-log(logSD)
    x.rawdensity<-dlnorm(Z,log=T)
    x.density<-x.rawdensity+scaling.log
    x.density[X==0]<--Inf
    lognorm.density.matrix<-matrix(exp(x.density),ncol=nrow(posterior.allchains),nrow=length(X))
  }
  density.quantile<-apply(lognorm.density.matrix,1,quantile,c(0.05,0.5,0.95))
  
  if(df==TRUE){
    density.quantile[,2]
    
    if(df){lognorm.density.df<-data.frame(duration=seq(0,range.eval,by=interval.eval),
                                          density=density.quantile[2,],
                                          density.95=density.quantile[3,],
                                          density.0.05=density.quantile[1,])
           
           return(lognorm.density.df)}
  }
    if(!df){return(rowMeans(lognorm.density.matrix))}
    
  }


# 
# sample.ndx<-c(250,2500,4750)
# density.quantile<-matrix(nrow=3,ncol=1001)
# Rprof(line.profiling=T)
# t1<-proc.time()
# for(i in 1:10){
# density.quantile<-apply(lognorm.density.matrix,1,quantile,c(0.05,0.5,0.95))
# #for(i in 1:nrow(lognorm.density.matrix)){density.quantile[,i]<-sort.int(lognorm.density.matrix[i,])[sample.ndx]}
# #tmp<-lognorm.plotdata.fun(btv.full.data$est,use.z=T)
# }
# t2<-proc.time()
# Rprof(NULL)
# t2-t1
# summaryRprof(lines="show")