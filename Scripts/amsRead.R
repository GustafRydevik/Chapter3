aws.sim.dir<-"/Users/gustafrydevik/simulations/Chapter3/aws/simulations"
library(ggplot2)
library(plyr)
library(data.table)
library(ggmcmc)
library(rjags)
library(coda)
# 
#  results.list<-vector(mode="list",length=length(dir(aws.sim.dir)))
# for(ndx in seq_along(dir((aws.sim.dir)))){
#   load.result<-tryCatch(load(dir(aws.sim.dir,full.names=T)[ndx]),error=function(e)e) ## should add a "try" here. 
#   if("error"%in%class(load.result)){next}
#   tmp<-ggs(as.mcmc.list(scenarios.results$est))
#   tmp.summary<-ddply(tmp,.(Parameter),summarize,mean=mean(value),
#                      lower.ci=quantile(value,0.025),
#                      upper.ci=quantile(value,0.975),
#                      filename=dir(aws.sim.dir)[ndx])
#   results.list[[ndx]]<-data.frame(subset(tmp.summary,Parameter%in%.(incidence,mean.incidence,sd,trend)),as.data.frame(scenarios.results$control.vars),
#                                   ndx=ndx)
#   print(ndx)
#   cat(paste0("\n",round(ndx/length(dir((aws.sim.dir)))*100),"%\n"))
# }
# 

#save(aws.results.list,file=file.path(output.path,"awsResList.RData"))
load(file=file.path(output.path,"awsResList.RData"))

 library(data.table)
 aws.diseaseType<-str_extract(dir(aws.sim.dir),"diseaseType[1:3]")
aws.ntest<-str_extract(dir(aws.sim.dir),"[a-zA-Z]+test")
aws.results.list<-lapply(1:length(aws.results.list),function(x){
   response<-tryCatch(data.frame(aws.results.list[[x]],diseaseType=aws.diseaseType[x],
                                 ntest=aws.ntest[x]),error=function(e)e)
   if("error"%in%class(response)){return(results.list[[x]])}else{
     return(response)
   }}
   )


aws.results.df<-rbindlist(aws.results.list)
library(reshape2)
aws.results.reshaped<-dcast(melt(aws.results.df,
                                 measure.vars=c("mean","lower.ci","upper.ci")
                                 ),
                            modelscript.name+n.infection+start.time+end.time+incidence+trend+Epi.model+diseaseType+ntest+ndx~Parameter+variable)
aws.results.reshaped<-(transform(aws.results.reshaped,ymax=trend_mean*end.time+incidence_mean))

aws.results.reshaped$row<-1:nrow(aws.results.reshaped)



aws.plot.data<-adply(aws.results.reshaped,1,function(x){
data.frame(x,est.line=c(x$incidence_mean,x$ymax),line.x=c(x$start.time,x$end.time),true.line=c(x$incidence,x$incidence+x$trend*x$end.time)
)})

##Add disease type, relative incidence

##Figure 3.9, but needs updating and subsetting!!!
library(RColorBrewer)
p.est.lines<-ggplot(data=subset(transform(aws.plot.data,rel.trend=trend/incidence),
                                          diseaseType=="diseaseType1"&rel.trend%in%c(0.05,0,-0.05)&incidence==0.01),
       aes(x=line.x,y=est.line,group=ndx,col=factor(ntest),group=ntest))+
  geom_abline(aes(intercept=incidence,slope=trend,col="True_line"),size=3,alpha=0.25)+
  geom_line()+
  xlab("Time")+ylab("Incidence")+
  ggtitle("Estimated trend lines \n by trend direction and sample size")+
  scale_colour_manual("Trend slope",values=c(True_line="black",Abtest=brewer.pal(3,"Set1")[1],
                                             NAtest=brewer.pal(3,"Set1")[2],
                                             Bothtest=brewer.pal(3,"Set1")[3]))+
  facet_grid(n.infection~rel.trend,scale="free")+theme_light(base_size=14)
ggsave_thesis(file.path(output.path,"fig_3_9_estimated_lines.png"),plot.name=p.est.lines)

##figure out a good xvalue here! 
ggplot(subset(transform(aws.results.reshaped,rel.trend=trend/incidence),
              diseaseType=="diseaseType1"&rel.trend%in%c(0.05,0,-0.05)),
       aes(x=ndx,y=trend_mean,
                                ymin=trend_lower.ci,
                                ymax=trend_upper.ci))+
  geom_pointrange()+geom_hline(y=0,col="red")+
  facet_grid(n.infection~rel.trend,scale="free")
  
ggplot(transform(aws.results.reshaped,rel.trend=trend/incidence),
       aes(x=ndx,y=trend_mean,
           ymin=trend_lower.ci,
           ymax=trend_upper.ci))+
  geom_pointrange()+geom_hline(y=0,col="red")+
  facet_grid(diseaseType+n.infection~rel.trend,scale="free")


##Plot of -0.05,0 and +0.05 trends, showing that trends are properly estimated with >1000 samples. 
##figure out a good ordering! 
aws.results.reshaped<-transform(aws.results.reshaped,rel.trend=trend/incidence)
aws.results.reshaped<-ddply(aws.results.reshaped,.(n.infection,rel.trend,incidence),function(df){data.frame(df,ndx.plot=1:nrow(df))})
ci.lines<-ggplot(subset(aws.results.reshaped,
              diseaseType=="diseaseType1"&rel.trend%in%c(0.05,0,-0.05)&incidence==0.01),
       aes(x=ntest,y=trend_mean,
           ymin=trend_lower.ci,
           ymax=trend_upper.ci),position="dodge")+
  geom_pointrange(size=1)+geom_hline(y=0,col="black")+
  facet_grid(n.infection~rel.trend)+
  xlab("Tests used")+ylab("Trend parameter value")+
  geom_point(aes(y=trend),col="red",size=3)+geom_point(aes(y=trend),col="white",size=2)+
  theme_light(base_size=14)+ggtitle("Posterior credible intervals for the trend\n by trend direction and sample size")
ggsave_thesis(file.path(output.path,"fig_3_10_posterior_ci_lines.png"),plot.name=ci.lines)

