
####Use 
#bash > nohup R --vanilla --slave --args seed 1000 burn.in 25000 adapt.iter 500 Set.endtime 63 nreps 200 < SerologyDataSim_Pertussis.R > nohup2.out 2>&1&
##to run in girion

###parameters

####Generating various simulated testing datasets, using functions from eventGenerators.R and Labdatagenerator.

##Parameters for using the batch package for multiple parallell runs


#Global parameters

got.mdrive<-length(dir("M:"))>0
is.win<-grepl("w32",R.Version()$platform)
is.girion<-sum(grep("girion",system("hostname",intern=TRUE)))
dropbox.sac<-"C:/Documents and Settings/GRydevik/My Documents/Dropbox"
dropbox.bioss<-"D:\\Dropbox"
dropbox.osx<-"/Users/gustafrydevik/Dropbox"
dropbox.girion<-"/home/gustaf/SeroNA_temp"
dropbox.path<-c(dropbox.osx,dropbox.sac,dropbox.bioss)[got.mdrive+is.win+1]
if(is.girion){dropbox.path<-dropbox.girion}
#sim.dir<-"D:/simulations/"
#sim.dir<-paste(dropbox.path,"/simulations",sep="")
sim.dir<-"/Users/gustafrydevik/simulations/Chapter3/"
if(is.girion)sim.dir<-file.path(dropbox.path,"simulations/")
##Project specific parameters
project.path<-file.path(dropbox.path,"PhD folder/Chapter3")
if(is.girion){project.path<-dropbox.girion}
data.path<-file.path(project.path,"Data")
script.path<-file.path(project.path,"Scripts")
output.path<-file.path(project.path,"Output")

lapply(dir(file.path(dropbox.path,"GusLib"),full.names=T),source)

##Loading linearLikelihood functions below here, which contains tsi.linear.li
lapply(dir(file.path(script.path,"hindcasting helper functions"),full.names=T),source)


##Note to self: pare this down substantially!
autolib(rjags)
autolib(batch)
autolib(reshape2)
autolib(ggplot2)
theme_thesis<-theme_minimal(base_size=14)
My.device<-Gen.device("png",res=400,width=12,height=3,units="in")


##Plotting the likelihood function here 
res<-100
a<-0.1
b<-0.8
c<-10
n.positive<-1000

  
incidence.range=c(0.1,0.01,0.001,0.00001)
trend.range<-rev(c(c(0.8,1/2,1/4,1/10),rev(-c(0.8,1/2,1/4,1/10)))) ##proportional final change at end of interval
list.ndx<-(match(a,incidence.range)-1)*length(trend.range)+match(b,trend.range)
b.abs<-a*b/c
prop.negative<-(1-(a+b.abs*c/2))^c
n.negative<-prop.negative/(1-prop.negative)*n.positive
plot.data<-na.omit(EndemicLinear(n.positive,
                                 incidence = a,
                                 trend = b.abs,start.time=0,end.time=c))

a.range<-seq(0.25*a,1.5*a,length.out=res)
b.range<-seq(max(-1/2*a/c,b.abs-abs(b.abs)/2),b.abs+2*abs(b.abs),length.out=res)
#b.range<-seq(-0.004,b.abs+3*abs(b.abs),length.out=res)


##This plots the value of tsi.linear.ll for values of a, and values of mean incidence
##Depending on the observed data in plot.data

ll.grid<-outer(a.range,b.range,function(X,Y){
  sapply(1:length(X),function(x){
    val<-tsi.linear.ll(pars =c(log(X[x]),log(1e-6+Y[x]+X[x]/c)),
                       data = plot.data,censorLimit. = c,popsize=n.positive+n.negative,
                       return.components=FALSE)
                       
    return(val)
  })
})
colnames(ll.grid)<-b.range
rownames(ll.grid)<-a.range

#ll.grid[ll.grid>20000]<-NA
#matplot(a.range,(ll.grid),type="l")

##Plot 1 in likelihood section
#Plot 3.3
##Add legend of true vs estimated maximum likelihood! 
tmp<-str_split(str_replace_all(paste(scan(file.path(data.path,"parula.txt"),what="char",sep = "\t"),collapse=""),"\\{|\\}|,|\\)","")," ")
tmp<-tmp[[1]][nzchar(tmp[[1]])]
parula_rgb<-paste("#",apply(matrix((
  as.character(as.hexmode(floor(as.numeric(tmp)*256)))),ncol=3,byrow = TRUE),
  1,paste0,collapse=""),sep="")


max.coord<-arrayInd(which.max(ll.grid),dim(ll.grid))
max.coord<-c(as.numeric(rownames(ll.grid)[max.coord[1]]),
             as.numeric(colnames(ll.grid)[max.coord[2]]))
p.ll.grid<-ggplot(transform(melt((ll.grid)),ll.quantile=ecdf(value)(value)),
                 aes(x=Var1,y=Var2))+
  geom_tile(aes(fill=(ll.quantile)))+
  scale_alpha_manual(values=c(0,1),guide="none")+
  #scale_fill_gradient(low="blue", high="yellow")+
  scale_fill_gradientn(colours=parula_rgb,"LL quantile")+
  geom_abline(slope=-0.2,intercept=0.028)+
  geom_point(x=a,y=b.abs,size=3,aes(col="True"))+
  geom_point(x=max.coord[1],y=max.coord[2],size=3,aes(col="Estimated"))+
  scale_colour_manual("LL maximum",values=c(Estimated="red",True="black"))+
  xlab(expression(alpha))+
  ylab(expression(beta))+#ylim(c(0.004,0.012))+
   #rownames(ll.grid)[arrayInd(which.min(ll.grid),dim(ll.grid))]
  ggtitle("Likelihood for different values of alpha and beta given data.
          \nLine indicates line of same average incidence")+theme_thesis
ggsave_thesis(file.path(output.path,"fig_3_1_loglik.png"),plot.name=p.ll.grid)
# 
# 
# 
# ##Zooming in on the same data 
# a.range2<-seq(a*0.7,a*1.3,length.out=res)
# b.range2<-seq(2*b.abs,-2*b.abs,length.out=res)
# ll.grid2<-outer(a.range2,b.range2,function(X,Y){
#   sapply(1:length(X),function(x){
#     tsi.linear.ll(pars =c(log(X[x]),log(1e-6+Y[x]+X[x]/c)),data = plot.data,censorLimit. = c,popsize =10000/(c*(a+b.abs*c/2)) )}
#   )})
# colnames(ll.grid2)<-b.range2
# rownames(ll.grid2)<-a.range2
# 
# ##Triangular for some reasons
# ggplot(melt((ll.grid2)),aes(x=Var1,y=Var2))+
#   geom_tile(aes(fill=value))+
#   scale_fill_gradient(low="red", high="green")+
#   geom_point(x=a,y=b.abs)
# 
# #ll.grid2[ll.grid2>-15000]<-NA
# 
# ggplot(melt((ll.grid2)),aes(x=Var1,y=Var2))+
#   geom_tile(aes(fill=value))+
#   scale_fill_gradient(low="red", high="green")+
#   geom_point(x=a,y=b.abs)+  #V2=0.004/0.12*V1; 
#   geom_point(x=a,y=tail(b.range2,1),col="blue")+
#   geom_point(x=a,y=head(b.range2,1),col="blue")+
#   geom_abline(slope=-0.2,intercept=0.028)  #Some sort of line of equivalence....
#                                            #Actually, this is just all lines passing through
#                                            #The mean incidence. suspect that none of the other parts differ much in that case.... 
# 
# hist(plot.data,freq=F,breaks=0:11,ylim=c(0,0.20))
# lines(0:11,lineardist(time = 0:11,alpha = a,beta=b.abs,chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = a,beta=b.abs,chain = 1)),col="red",lwd=2)
# for(atmp in sample(a.range2,10)){
#   for(btmp in sample(b.range2,10)){
#     lines(0:11,lineardist(time = 0:11,alpha = atmp,beta=btmp,chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = atmp,beta=btmp,chain = 1)),
#           lwd=0.25,col=adjustcolor("black",alpha.f=0.5))
#   }
# }
# lines(0:11,lineardist(time = 0:11,alpha = a,beta=tail(b.range2,1),chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = a,beta=tail(b.range2,1),chain = 1)),col="blue")
# lines(0:11,lineardist(time = 0:11,alpha = a,beta=head(b.range2,1),chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = a,beta=head(b.range2,1),chain = 1)),col="blue")
# 
# 
# 
# 
# 
# incidence.range=c(0.1,0.01,0.001,0.00001)
# trend.range<-rev(c(c(0.8,1/2,1/4,1/10),rev(-c(0.8,1/2,1/4,1/10)))) ##proportional final change at end of interval
# names(incidence.range)<-c(0.1,0.01,0.001,0.00001)
# names(trend.range)<-trend.range
# 
# mean.incidence.range=c(0.1,0.01,0.001,0.00001)
# 
# ll.image.list<-vector(mode="list",length=length(incidence.range)*length(trend.range))
# for(a in incidence.range[1]){
#   for(b in trend.range[1]){
#     list.ndx<-(match(a,incidence.range)-1)*length(trend.range)+match(b,trend.range)
#     b.abs<-a*b/c
# 
#     
#     list.ndx<-(match(a,incidence.range)-1)*length(trend.range)+match(b,trend.range)
#     b.abs<-a*b/c
#     plot.data<-na.omit(EndemicLinear(1000/(c*(a+b.abs*c/2)),
#                                      incidence = a,
#                                      trend = b.abs,start.time=0,end.time=c))
#     
#     a.range<-seq(1/2*a,1.5*a,length.out=res)
#     b.range<-seq(max(-1/2*a/c,b.abs-abs(b.abs)/2),min(1/2*a/c,b.abs+abs(b.abs)/2),length.out=res)
#     ll.grid<-outer(a.range,b.range,function(X,Y){
#       sapply(1:length(X),function(x){
#         tsi.linear.ll(pars =c(log(X[x]),log(1e-6+Y[x]+X[x]/c)),data = plot.data,censorLimit. = c,
#                       popsize =1000/(c*(a+b.abs*c/2)) )}
#       )
#       colnames(ll.grid)<-b.range
#       rownames(ll.grid)<-a.range
#       ll.image.list[[list.ndx]]<-ll.grid  
#     })
#   }
# }
# 
# ## Evaluating performance here
# 
# 
# 
# 
# trendbias.grid<-outer(incidence.range,trend.range,function(X,Y){
#   sapply(1:length(X),function(x){
#     tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c,n.obs=3000)$trend.bias}
#     )
#   })
# 
# trendest.grid<-outer(incidence.range,trend.range,function(X,Y){
#     sapply(1:length(X),function(x)tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c)$trend.est)
#   })
# 
# incidence.est.grid<-outer(incidence.range,trend.range,function(X,Y){
#   sapply(1:length(X),function(x)tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c)$incidence.est)
# }
# )
# 
# incidence.bias.grid<-outer(incidence.range,trend.range,function(X,Y){
#   sapply(1:length(X),function(x)tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c)$incidence.bias)
# }
# )
# est.df<-merge(melt(incidence.est.grid),melt(trendest.grid),by=c("Var1","Var2"))
# names(est.df)<-c("true.inc","rel.trend","est.inc","est.trend")
# est.df$true.trend<-est.df$true.inc*est.df$rel.trend/c
# est.df$trend.relbias<-est.df$est.trend/est.df$true.trend
# matplot(trend.range,t(trendbias.grid),type="l")
# 
# matplot(trend.range,t(trendest.grid),type="l")
# 
# melt(trendest.grid)
# ggplot(melt(incidence.bias.grid),aes(y=value,x=Var2,group=Var1,col=factor(Var1)))+geom_line()
# ggplot(melt(trendest.grid),aes(y=value/(Var2*Var1/c),x=Var2))+geom_line()+facet_wrap(~Var1,scales="free")
# ggplot(melt(trendest.grid),aes(y=value,x=Var2,group=Var1,col=factor(Var1)))+geom_line()
# 
# 
# 
# 
# ##tsi match
# estTrend.df<- as.data.frame(apply(
#   est.df,1,function(x)lineardist(0:10,x[3],x[4])))
# colnames(estTrend.df)<-gsub("-","minus",paste(format(est.df$true.inc,scientific=FALSE,drop0trailing=TRUE),est.df$rel.trend,sep=".."))
# estTrend.df<-data.frame(Time=0:10,estTrend.df)
# estTrend.df<-melt(estTrend.df,id.vars="Time",variable.name="combo")
# estTrend.df$true.inc<-as.numeric(gsub("X","",sapply(strsplit(as.character(estTrend.df$combo),"\\.\\."),"[",1)))
# estTrend.df$rel.trend<-as.numeric(gsub("minus","-",sapply(strsplit(as.character(estTrend.df$combo),"\\.{2,}"),"[",2)))
# 
# ggplot(estTrend.df,aes(x=Time,y=value))+geom_line()+facet_grid(true.inc~rel.trend)+
#   geom_histogram(data=rlinear.tsi(1000,incidence=a,trend=b),aes(y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=1)
# 
# 
# ###linear match
# 
# estTrend.df<- as.data.frame(apply(
#   est.df,1,function(x)x[3]+x[4]*(0:10)))
# colnames(estTrend.df)<-gsub("-","minus",paste(format(est.df$true.inc,scientific=FALSE,drop0trailing=TRUE),est.df$rel.trend,sep=".."))
# estTrend.df<-data.frame(Time=0:10,estTrend.df)
# estTrend.df<-melt(estTrend.df,id.vars="Time",variable.name="combo")
# estTrend.df$true.inc<-as.numeric(gsub("X","",sapply(strsplit(as.character(estTrend.df$combo),"\\.\\."),"[",1)))
# estTrend.df$rel.trend<-as.numeric(gsub("minus","-",sapply(strsplit(as.character(estTrend.df$combo),"\\.{2,}"),"[",2)))
# 
# ggplot(estTrend.df,aes(x=Time,y=value))+geom_line()+facet_grid(true.inc~rel.trend)+geom_abline(aes(intercept=true.inc,slope=true.inc*rel.trend/10),col="red")
# 
#                                   
# 
# ##Looking at plotting the likelihood function 
# 
# #This seems to indicate that the best estimate is biased and negative... 
# test.data<-EndemicLinear(n.infections =10000,incidence = 0.01,trend=0.01/(30*2),start.time = 0,end.time = 30)
# plot(seq(-0.01/30,0.01/30,length.out=100),
#      sapply(seq(-0.01/30,0.01/30,length.out=100),
#             function(x)tsi.linear.ll(c(log(0.01),log(x+2*0.01/30)),data=test.data,popsize=10000,censorLimit. = 30)
#             )
#      )
# 
# min.ndx<-which.min(sapply(seq(-0.01/30,0.01/30,length.out=100),
#                  function(x)tsi.linear.ll(c(log(0.01),log(x+2*0.01/30)),data=test.data,popsize=10000,censorLimit. = 30)))
# abline(v=0.01/(30*2),col="red")
# abline(v=seq(-0.01/30,0.01/30,length.out=100)[min.ndx],col="blue")
# seq(-0.01/30,0.01/30,length.out=100)[min.ndx]
# 
# hist(test.data,freq=FALSE,ylim=c(0,0.1))
# for(beta in seq(-0.01/30,0.01/30,length.out=10)){lines(dlinear.tsi(InfTime =0:30,censorLimit = 30,incidence = 0.01,trend = beta),col="black",lwd=0.5)}
# lines(dlinear.tsi(InfTime =0:30,censorLimit = 30,incidence = 0.01,trend = 0.005/30),col="red")
# lines(dlinear.tsi(InfTime =0:30,censorLimit = 30,incidence = 0.01,trend = seq(-0.01/30,0.01/30,length.out=100)[min.ndx]),col="blue")
# ##Biased downards because there is no individuals infected at T==30...
