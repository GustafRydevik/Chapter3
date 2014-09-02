
hindcast.performance.eval<-function(sim.name,ma.width=7,scenario.pars=NULL,...){
  require(stringr,rjags)
  if(is.character(sim.name)){simresults.name<-load(sim.name)
  results.ndx<-which(sapply(lapply(scenarios.results,names),function(x)any(grepl("est",x))))
  epi.obs<-as.numeric(str_extract(str_extract(sim.name,"endtime\\_[[:digit:]]+"),"[[:digit:]]+"))
  epi.start<-as.numeric(str_extract(str_extract(sim.name,"start|starttime\\_[[:digit:]]+"),"[[:digit:]]+"))
  pert.or.btv<-str_extract(tail(strsplit(sim.name,"/")[[1]],1),"(pert)|(btv)")
  ntests<-str_extract(sim.name,"[[:alpha:]]+test")
  sample.size<-nrow(get(simresults.name)[[results.ndx]]$true.values)
  #sample.size<-str_extract(str_extract(sim.name,"samplesize[[:digit:]]+"),"[[:digit:]]+")
  sim.results<-get(simresults.name)[[results.ndx]]$est
  }else{
    sim.results<-sim.name
    for(i in names(scenario.pars)){
      assign(i,scenario.pars[[i]])
    }
  }
  
  if(pert.or.btv=="pert"|pert.or.btv=="pertussis"){true.curve<-pertussis.wi.Fun(1,start.time=epi.start,end.time=epi.obs,Actual=T)}
  if(pert.or.btv=="btv"){true.curve<-btv.2008.Fun(1,start.time=epi.start,end.time=epi.obs,Actual=T)}
  
  true.density<-hist(true.curve,breaks=seq(0,epi.obs-epi.start,by=1),plot=F)$density
  ###Removing the NA's in the below.
  true.density.smooth<-c(cumsum(true.density[1:((ma.width-1)/2)])/(1:((ma.width-1)/2)),
                         na.omit(filter(true.density,rep(1,ma.width)/ma.width,sides=2)),
                         (cumsum(tail(true.density,ma.width))/(1:(ma.width)))[((ma.width-1)/2+2):ma.width])
  estimated.density<-lognorm.plotdata.fun(sim.results,df=FALSE,interval.eval = 1,range.eval=(epi.obs-epi.start)*4)
  estimated.density.df<-lognorm.plotdata.fun(sim.results,df=TRUE,interval.eval = 1,range.eval=(epi.obs-epi.start)*4)
  
  est.r2<-cor((true.density),estimated.density[2:(epi.obs-epi.start+1)])^2
  est.r2.smooth<-cor(true.density.smooth,estimated.density[2:(epi.obs-epi.start+1)])^2
  #est.epistart<-do.call("c",sim.results[,"EpiStart",])
  est.peaktime<-do.call("c",sim.results[,"peak.time",])
  est.duration<-do.call("c",sim.results[,"duration",])
  est.epistart<-qlnorm(1-0.5/sum(true.curve),log(est.peaktime),log(est.duration))
  peaktime.50<-quantile(est.peaktime,0.5)
  peaktime.05<-quantile(est.peaktime,0.05)
  peaktime.95<-quantile(est.peaktime,0.95)
  duration.50<-quantile(est.duration,0.5)
  duration.05<-quantile(est.duration,0.05)
  duration.95<-quantile(est.duration,0.95)  
  epistart.50<-quantile(est.epistart,0.5)
  epistart.95<-quantile(est.epistart,0.95)
  epistart.05<-quantile(est.epistart,0.5)
  epistart.bias.rel<-quantile((est.epistart)/(epi.obs-epi.start),c(0.05,0.5,0.95))[2]
  epistart.bias<-quantile(abs((est.epistart)-(epi.obs-epi.start)),c(0.05,0.5,0.95)) ## change this to non-absolute...
  epistart.bias.50<-epistart.bias[2]
  epistart.bias.05<-epistart.bias[1]
  epistart.bias.95<-epistart.bias[3]
  
  est.absmean<-(mean(abs(sample.size*estimated.density[2:(epi.obs-epi.start+1)]-sample.size*true.density.smooth)))
  est.rmsep<-sqrt(mean((sample.size*estimated.density[2:(epi.obs-epi.start+1)]-sample.size*true.density.smooth)^2))
  #gelman.convergence<-gelman.diag(sim.results[,!(str_detect(colnames(sim.results[[1]]),"(InfTime)")),])
  gelman.convergence<-gelman.diag(sim.results[,c("peak.time","duration"),])
  results.df<-data.frame(est.r2=est.r2,
                         est.r2.smooth=est.r2.smooth,
                         est.rmsep=est.rmsep,
                         epistart.50=epistart.50,
                         epistart.05=epistart.05,
                         epistart.95=epistart.95,
                         epistart.bias.50=epistart.bias.50,
                         epistart.bias.05=epistart.bias.05,
                         epistart.bias.95=epistart.bias.95,
                         epistart.relative.bias=epistart.bias.rel,
                         duration.50=duration.50,
                         duration.05=duration.05,
                         duration.95=duration.95,
                         peaktime.50=peaktime.50,
                         peaktime.05=peaktime.05,
                         peaktime.95=peaktime.95,
                         overall.convergence=gelman.convergence$mpsrf,
                         epi.obs=epi.obs,
                         epi.start=epi.start,
                         pathogen=pert.or.btv,one.or.two=ntests,
                         sample.size=sample.size)
  return(list(results.df=results.df,
              detailed.results=list(
                est.r2=est.r2,bias=epistart.bias,
                gelman.convergence=gelman.convergence,
                density=estimated.density.df)))
}
