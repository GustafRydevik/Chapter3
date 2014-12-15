library(ggplot2)
#### example of a Lotka-Volterra graph 
lv.example.df<-data.frame(time=seq(1,10,by=0.1),LotkaVolterra.Fun(diseaseType1)(seq(1,10,by=0.1)))
lv.example.p<-ggplot(data=lv.example.df,aes(x=time))+
  geom_line(aes(y=t1),col="red",size=2)+
  geom_line(aes(y=t2),col="blue",size=2)+
  ylab("test response")+ggtitle("Example of a Lotka-Volterra modelled test response curve")+theme_minimal()
ggsave(file.path(output.path,"lV_example.png"),plot=lv.example.p,width=8,height=5,units="in",dpi=round(400/7))


lv.pars.df<-data.frame(
  Type=c("Incubating pathogen","Fast-acting pathogen","Chronic infection with acute phase"),
  b_NA=c(diseaseType1$parms["growth.na"],diseaseType2$parms["growth.na"],diseaseType3$parms["growth.na"]),
  h=c(diseaseType1$parms["killed.na"],diseaseType2$parms["killed.na"],diseaseType3$parms["killed.na"]),
  b_ab=c(diseaseType1$parms["feed.ab"],diseaseType2$parms["feed.ab"],diseaseType3$parms["feed.ab"]),
  d_ab=c(diseaseType1$parms["dieoff.ab"],diseaseType2$parms["dieoff.ab"],diseaseType3$parms["dieoff.ab"])
)
library(pander)
pandoc.table(lv.pars.df)

###Lotka Volterra graphs for three types of diseases. 
lv.all.df<-data.frame(Type=rep(c("Incubating pathogen",
                                 "Fast-acting pathogen","Chronic infection with acute phase"),
                               each=length(seq(1,20,by=0.01))),
                      time=seq(1,20,by=0.01),ab=c(
                        LotkaVolterra.Fun(diseaseType1)(seq(1,20,by=0.01))[,1],
                        LotkaVolterra.Fun(diseaseType2)(seq(1,20,by=0.01))[,1],
                        LotkaVolterra.Fun(diseaseType3)(seq(1,20,by=0.01))[,1]),
                      na=c(
                        LotkaVolterra.Fun(diseaseType1)(seq(1,20,by=0.01))[,2],
                        LotkaVolterra.Fun(diseaseType2)(seq(1,20,by=0.01))[,2],
                        LotkaVolterra.Fun(diseaseType3)(seq(1,20,by=0.01))[,2])
)


lv.timeplot.all<-ggplot(data=lv.all.df,aes(x=time))+
  geom_line(aes(y=ab),col="red",size=2)+
  geom_line(aes(y=na),col="blue",size=2)+facet_wrap(~Type,ncol=1,scales="free_y",as.table=FALSE)+
  scale_x_continuous(limits=c(0,20))+ylab("test value")

lv.phaseplot.all<-ggplot(data=lv.all.df,aes(x=na,y=ab))+
  geom_point(col="red",size=2)+
  facet_wrap(~Type,as.table=FALSE,scales="free",ncol=1)+
  ylab("Nucleic acid")+xlab("Antibodies")
lv.phaseplot.all
library(gridExtra)
lv.time.phase<-arrangeGrob(lv.timeplot.all,lv.phaseplot.all,nrow=1)
ggsave(file.path(output.path,"LV_phase_time.png"),plot=lv.time.phase,width=8,height=8,units="in",dpi=round(400/7))



##Posterior plots
library(ggmcmc)


##Evaluating likelihood

hindcast.ll<-function(data,type="linear",diseaseType="DiseaseType1",pars){
  if(type!=linear){error("Only the linear model is implemented!")}
  if(type==linear){
    trend<-pars$trend
    incidence<-pars$incidence
    sd<-pars[grep("sd",names(pars))]
    InfTimes<-pars[grep("InfTime",names(pars))]
    
    E.data<-LotkaVolterra.Fun(disease=get(diseaseType))(InfTimes)
    data.diff<-data/E.data
    plnorm()
  }
  
}



load(file.path(sim.dir,"diseaseType3_EndemicDecreaseFixedSD_seed_1002_t2_samplesize5000_00_gr1.RData"))
#load(file.path(sim.dir,"diseaseType3_EndemicKnownTimes_seed_1002_t2_samplesize2000_00_gr1.RData"))

censorLimit<-(scenarios.results$control.vars$end.time-scenarios.results$control.vars$start.time)
ggs.decrease<-ggs(scenarios.results$est)
decrease.true.inftime<-data.frame(true.inftime=scenarios.results$true.values$infection.times)
decrease.true.inftime$Parameter<-paste0("InfTime[",1:nrow(decrease.true.inftime),"]")
decrease.true.inftime$scenario<-"decrease"
decrease.true.pars<-data.frame(scenarios.results$control.vars,scenario="decrease")



load(file.path(sim.dir,"diseaseType3_EndemicConstantFixedSD_seed_1000_t2_samplesize5000_00_gr1.RData"))
#load(file.path(sim.dir,"diseaseType3_EndemicKnownTimes_seed_1000_t2_samplesize2000_00_gr1.RData"))
ggs.constant<-ggs(scenarios.results$est)

constant.true.inftime<-data.frame(true.inftime=scenarios.results$true.values$infection.times)
constant.true.inftime$Parameter<-paste0("InfTime[",1:nrow(constant.true.inftime),"]")
constant.true.inftime$scenario<-"constant"
constant.true.pars<-data.frame(scenarios.results$control.vars,scenario="constant")


load(file.path(sim.dir,"diseaseType3_EndemicIncreaseFixedSD_seed_1001_t2_samplesize5000_00_gr1.RData"))
#load(file.path(sim.dir,"diseaseType3_EndemicKnownTimes_seed_1001_t2_samplesize2000_00_gr1.RData"))
ggs.increase<-ggs(scenarios.results$est)
increase.true.inftime<-data.frame(true.inftime=scenarios.results$true.values$infection.times)
increase.true.inftime$Parameter<-paste0("InfTime[",1:nrow(increase.true.inftime),"]")
increase.true.inftime$scenario<-"increase"
increase.true.pars<-data.frame(scenarios.results$control.vars,scenario="increase")


ggs.all<-rbind(data.frame(ggs.increase,scenario="increase"),
               data.frame(ggs.decrease,scenario="decrease"),
               data.frame(ggs.constant,scenario="constant"))

true.pars<-rbind(increase.true.pars,decrease.true.pars,constant.true.pars)
true.pars<-transform(true.pars,mean.incidence=incidence+end.time/2*trend)
est.trendpars<-subset(ggs.all,Parameter%in%c("incidence","trend","mean.incidence"))
est.inftime<-subset(ggs.all,grepl("InfTime",Parameter))

trendpars.mean<-ddply(subset(est.trendpars,max(Iteration)-Iteration<101),.(Parameter,Chain,scenario),summarise,meanPar=mean(value))
#true.inftime<-data.frame(true.inftime=scenarios.results$true.values$infection.times)
#true.inftime$Parameter<-paste0("InfTime[",1:nrow(true.inftime),"]")


library(data.table)
increase.true.inftime<-data.table(increase.true.inftime, key="Parameter,scenario")
decrease.true.inftime<-data.table(decrease.true.inftime, key="Parameter,scenario")
constant.true.inftime<-data.table(constant.true.inftime, key="Parameter,scenario")
est.inftime<-data.table(est.inftime,key="Parameter,scenario")

all.inftime<-merge(rbindlist(list(increase.true.inftime,
                         decrease.true.inftime,
                         constant.true.inftime)),est.inftime,by=c("Parameter","scenario"))
#mean.inftime<-data.frame(true.inftime,ddply(ggs.inftime,.(Parameter),summarise,meanInfTime=mean(value)))


ggplot(data=est.trendpars,
       aes(x=value,group=Chain))+
  geom_density(aes(fill=factor(Chain)),alpha=0.3)+
  facet_grid(scenario~Parameter,scales="free")+
  geom_vline(data=subset(melt(true.pars,variable.name="Parameter"),Parameter%in%c("incidence","trend","mean.incidence")),
             aes(xintercept=value))


#Residual plot below - indicates that for InfTime>17, the probability that the time of infection is ~0
# is non-zero


lineardist<-function(time,alpha,beta,chain){
  alpha<-alpha[chain]
  beta<-beta[chain]
  (alpha+beta*time)*exp(-(alpha+beta*time/2)*time)}

estTrend.df<- as.data.frame(apply(
                          recast(trendpars.mean,Chain+scenario~Parameter,measure.var="meanPar"),
                          1,function(x)lineardist(1:censorLimit,as.numeric(x[3]),as.numeric(x[5]),1)))

trueTrend.df<-do.call("rbind",
        apply(true.pars,1,function(x){
          data.frame(scenario=x[8],Time=1:censorLimit,value=lineardist(1:censorLimit,as.numeric(x[5]),as.numeric(x[6]),1))}
          ))
trueTrend.df<-ddply(trueTrend.df, .(scenario), mutate,
                   weight = sum(value))
colnames(estTrend.df)<-paste0("Chain",rep(1:5,each=3),".scenario.",rep(c("increase","decrease","constant"),5))
estTrend.df<-data.frame(Time=1:censorLimit,estTrend.df)
estTrend.df<-melt(estTrend.df,id.vars="Time",variable.name="Chain")
estTrend.df$scenario<-sapply(strsplit(as.character(estTrend.df$Chain),"\\."),"[",3)
estTrend.df$Chain<-sapply(strsplit(as.character(estTrend.df$Chain),"\\."),"[",1)

estTrend.df$Chain<-as.numeric(as.factor(estTrend.df$Chain))
estTrend.df$TrueCurve<-lineardist(estTrend.df$Time,1/100,1/100/50,1)
estTrend.df<-ddply(estTrend.df, .(Chain,scenario), mutate,
      weight = sum(value))

ggplot(all.inftime,aes(x=factor(round(true.inftime)),y=value-true.inftime))+
  geom_violin()+
  facet_grid(scenario~Chain)

ggplot(all.inftime,aes(x=value))+geom_histogram(aes(y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=1)+
  facet_grid(scenario~Chain)+
  geom_line(data=estTrend.df,aes(x=Time,y=value/weight),size=1.5)+
  geom_line(data=trueTrend.df,aes(x=Time,y=value/weight),col="red",alpha=0.5,size=1.5)
  
ggplot(all.inftime,aes(x=true.inftime))+geom_histogram(aes(y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=1)+
  facet_grid(scenario~Chain)+
  geom_line(data=estTrend.df,aes(x=Time,y=value/weight),size=1.5)+
  geom_line(data=trueTrend.df,aes(x=Time,y=value/weight),col="red",alpha=0.5,size=1.5)


###scatterplot of trend vs incidence
ggplot(recast(est.trendpars,scenario+Chain+Iteration~Parameter,measure.var="value"),aes(x=incidence,y=trend,col=Iteration))+
  geom_line(alpha=0.75)+geom_point(alpha=0.75)+
  geom_vline(data=true.pars,
             aes(xintercept=incidence),alpha=0.5,size=2)+
  geom_hline(data=true.pars,
                           aes(yintercept=trend),alpha=0.5,size=2)+
             facet_grid(scenario~Chain,scales="free") 



##scatterplot for evaluating convergence
ggplot(transform(recast(est.trendpars,scenario+Chain+Iteration~Parameter,measure.var="value"),
                 mean.incidence=incidence+trend*censorLimit/2),aes(x=mean.incidence,y=trend))+
  geom_vline(data=true.pars,
             aes(xintercept=trend*censorLimit/2+incidence),alpha=0.5,size=2)+
  geom_hline(data=true.pars,
             aes(yintercept=trend),alpha=0.5,size=2)+
  geom_point(alpha=0.4,aes(col=factor(Chain)))+#geom_line(alpha=0.75)+
  geom_point(data=true.pars,
             aes(y=trend,x=trend*censorLimit/2+incidence,col=scenario),size=5)+
  facet_grid(scenario~Chain,scales="free") 




## Scatterplot of the posterior for results purposes 
ggplot(transform(recast(est.trendpars,scenario+Chain+Iteration~Parameter,measure.var="value"),
   mean.incidence=incidence+trend*censorLimit/2),aes(x=mean.incidence,y=trend))+
  geom_vline(data=true.pars,
             aes(xintercept=trend*censorLimit/2+incidence),alpha=0.5,size=2)+
  geom_hline(data=true.pars,
             aes(yintercept=trend),alpha=0.5,size=2)+
   geom_point(alpha=0.4,aes(col=scenario))+#geom_line(alpha=0.75)+
  geom_point(data=true.pars,
              aes(y=trend,x=trend*censorLimit/2+incidence,col=scenario),size=5)


##True histogram 
#hist(true.inftime$true.inftime,freq=FALSE)
#lines(1:30,lineardist(1:30,0.01,0.01/50,chain=1)/mean(1-is.na(true.inftime$true.inftime)))
