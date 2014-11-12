
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



