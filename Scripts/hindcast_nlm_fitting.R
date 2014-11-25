
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
lapply(dir(file.path(script.path,"multitest library"),full.names=T),source)

##Note to self: pare this down substantially!
autolib(rjags)
autolib(batch)
autolib(reshape2)

My.device<-Gen.device("png",res=400,width=12,height=3,units="in")

#### Reading in  functions specific for this project 


lineardist<-function(time,alpha,beta,chain){
  alpha<-alpha[chain]
  beta<-beta[chain]
  (alpha+beta*time)*exp(-(alpha+beta*time/2)*time)}


dlinear.tsi<-function(InfTime,censorLimit=30,incidence=0.01,trend=0.01/(30*2)){
  area<-(1-exp(-(incidence*censorLimit+trend*censorLimit^2/2)))
  
  pI<-(incidence+trend*InfTime)*exp(-(incidence+trend*(InfTime)/2)*InfTime)/area*(censorLimit>InfTime)
  return(pI)
}


### cumulative probablity of linear.tsi 
rlinear.tsi<-function(n=1000,censorLimit=30,incidence=0.1,trend=0.1/(30*2)){
  if(trend<(-incidence/censorLimit)){stop(simpleError("Too steep linear trend results in \nnegative incidence over the interval"))}
PI.lookup<-data.frame(I=seq(0,censorLimit,by=0.1),PI=sapply(seq(0,censorLimit,by=0.1),
                                                 function(x){
                                                   integrate(dlinear.tsi,0,x,incidence=incidence,trend=trend)$value}
                                                 ))
require(data.table)
runif.dt<-data.table(data.frame(P=runif(n=n)),key="P")
result.linear.tsi<-data.table(PI.lookup,key="PI")[runif.dt,list(P,I),roll="nearest"]
return(result.linear.tsi$I)
}

a=0.01
censorLimit=30
b=0.01/30

###Why is the variability much higher in the "direct" generator, than in the simulated one 
par(mfrow=c(2,1))
hist(round(rlinear.tsi(1000,incidence=a,trend=b)),breaks=0:31)
hist(EndemicLinear(round(1000/(c*(a+b*c/2))),incidence=a,trend=b,start.time=0,end.time=30),breaks=0:31)

plot(sapply(seq(0,30,by=0.1),function(x)integrate(dlinear.tsi,x,30)$value))
test.data1<-rlinear.tsi(1000,incidence=a,trend=b,censorLimit=censorLimit)
test.data2<-EndemicLinear(round(1000/(censorLimit*(a+b*censorLimit/2))),incidence=a,trend=-b,start.time=0,end.time=censorLimit)


#A function for estimating the RMSEP between a the linear trend curve decided by pars, and the histogram of data
tsi.linear.rmse<-function(pars,censorLimit.=censorLimit,data=obs.data,transform.pars=TRUE){
  incidence.<-pars[1]
  trend.<-pars[2]
  if(transform.pars){
    incidence.<-exp(incidence.)
    trend.<-exp(trend.)-incidence./censorLimit.
  }
  data.prop<-hist(data,breaks=0:censorLimit.,plot=F)
  predprop<-dlinear.tsi(data.prop$mids,incidence=incidence.,trend=trend.,censorLimit=censorLimit.)
  predprop.se<-(predprop-data.prop$density)^2
  return(sqrt(mean(predprop.se)))
}

#A function for evaluation the data likelihood of the linear trend curve decided by pars, given data.
#Accepts parameters>0, ie logtransformed, so as to work well in conjunction with nlm

tsi.linear.ll<-function(pars,censorLimit.,data=obs.data,transform.pars=TRUE,popsize){
  incidence.<-pars[1]
  trend.<-pars[2]
  if(transform.pars){
    incidence.<-exp(incidence.)
    trend.<-exp(trend.)-incidence./censorLimit.
  }
  InfTime<-data[!is.na(data)]
  mean.incidence<-incidence.+trend.*censorLimit./2
  outside.range<-any(incidence.-trend.*InfTime<0)
  trend.prior.ll<-dunif(trend.,-2*mean.incidence/censorLimit.,2*mean.incidence/censorLimit.,log=T)
  mean.inc.ll<-dbeta(1-(1-mean.incidence)^censorLimit.,length(InfTime)/censorLimit.+1,(popsize-length(InfTime))/censorLimit.+1,log=T)
  area<-log((1-exp(-(incidence.+trend.*censorLimit./2)*censorLimit.)))
  loglik<-trend.prior.ll+mean.inc.ll-area+
  log((incidence.+trend.*InfTime)*exp(-(incidence.+trend.*(InfTime)/2)*InfTime))+log((1-outside.range))
  return(-sum(loglik))
}

#A function for attempting to use nlm to fit the linear tsi curve. 
tsi.linear.est<-function(incidence=0.01,
                         trend=0.01/30,
                         censorLimit=30,
                         simtype="population",
                         n.obs=1000,engine="nlm",transform.pars=TRUE,type="ll"){
  
  
  approx.popsize<-round(n.obs/(censorLimit*(incidence+trend*censorLimit/2)))
  obs.data<-if(simtype=="population"){EndemicLinear(approx.popsize,incidence=incidence,trend=trend,start.time=0,end.time=censorLimit)
                                 }else{rlinear.tsi(n.obs,incidence=incidence,trend=trend,censorLimit=censorLimit)}
 
  obs.data<-obs.data[!is.na(obs.data)]
  prior.inc<-exp(censorLimit*log(1-length(obs.data)/(approx.popsize)))
  pars.est<-nlm(tsi.linear.ll,p=c(log(prior.inc),
                                  log((prior.inc)/(censorLimit.))),
                transform.pars=transform.pars,data=obs.data,popsize=approx.popsize,censorLimit.=censorLimit,
                steptol=1e-8)
  if(transform.pars){
  pars.est$estimate[1]<-exp(pars.est$estimate[1])
  pars.est$estimate[2]<-exp(pars.est$estimate[2])-pars.est$estimate[1]/censorLimit
  }
  return(list(incidence.est=pars.est$estimate[1],
              trend.est=pars.est$estimate[2],
              incidence.bias=incidence-pars.est$estimate[1],
              trend.bias=trend-pars.est$estimate[2])
  )
}





##Plotting the likelihood function here 
res<-100
a<-0.1
b<-0.8
c<-10
list.ndx<-(match(a,incidence.range)-1)*length(trend.range)+match(b,trend.range)
b.abs<-a*b/c
plot.data<-na.omit(EndemicLinear(10000/(c*(a+b.abs*c/2)),
                                 incidence = a,
                                 trend = b.abs,start.time=0,end.time=c))

a.range<-seq(1/4*a,2*a,length.out=res)
b.range<-seq(max(-1/2*a/c,b.abs-abs(b.abs)/2),b.abs+2*abs(b.abs),length.out=res)
ll.grid<-outer(a.range,b.range,function(X,Y){
  sapply(1:length(X),function(x){
    tsi.linear.ll(pars =c(log(X[x]),log(1e-6+Y[x]+X[x]/c)),data = plot.data,censorLimit. = c,popsize =10000/(c*(a+b.abs*c/2)) )}
  )})
colnames(ll.grid)<-b.range
rownames(ll.grid)<-a.range

ll.grid[ll.grid>20000]<-NA

ggplot(melt((ll.grid)),aes(x=Var1,y=Var2))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low="red", high="green")+
  geom_point(x=a,y=b.abs)

##Zooming in on the same data 
a.range2<-seq(a*0.7,a*1.5,length.out=res)
b.range2<-seq(2*b.abs,-2*b.abs,length.out=res)
ll.grid2<-outer(a.range2,b.range2,function(X,Y){
  sapply(1:length(X),function(x){
    tsi.linear.ll(pars =c(log(X[x]),log(1e-6+Y[x]+X[x]/c)),data = plot.data,censorLimit. = c,popsize =10000/(c*(a+b.abs*c/2)) )}
  )})
colnames(ll.grid2)<-b.range2
rownames(ll.grid2)<-a.range2


ggplot(melt((ll.grid2)),aes(x=Var1,y=Var2))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low="red", high="green")+
  geom_point(x=a,y=b.abs)

ll.grid2[ll.grid2>-15000]<-NA

ggplot(melt((ll.grid2)),aes(x=Var1,y=Var2))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient(low="red", high="green")+
  geom_point(x=a,y=b.abs)+  #V2=0.004/0.12*V1; 
  geom_point(x=a,y=tail(b.range2,1),col="blue")+
  geom_point(x=a,y=head(b.range2,1),col="blue")+
  geom_abline(slope=-0.2,intercept=0.028)  #Some sort of line of equivalence....
                                           #Actually, this is just all lines passing through
                                           #The mean incidence. suspect that none of the other parts differ much in that case.... 

hist(plot.data,freq=F,breaks=0:11,ylim=c(0,0.20))
lines(0:11,lineardist(time = 0:11,alpha = a,beta=b.abs,chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = a,beta=b.abs,chain = 1)),col="red",lwd=2)
for(atmp in sample(a.range2,10)){
  for(btmp in sample(b.range2,10)){
    lines(0:11,lineardist(time = 0:11,alpha = atmp,beta=btmp,chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = atmp,beta=btmp,chain = 1)),
          lwd=0.25,col=adjustcolor("black",alpha.f=0.5))
  }
}
lines(0:11,lineardist(time = 0:11,alpha = a,beta=tail(b.range2,1),chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = a,beta=tail(b.range2,1),chain = 1)),col="blue")
lines(0:11,lineardist(time = 0:11,alpha = a,beta=head(b.range2,1),chain = 1)*10/sum(lineardist(time = seq(0,10,by=0.1),alpha = a,beta=head(b.range2,1),chain = 1)),col="blue")





incidence.range=c(0.1,0.01,0.001,0.00001)
trend.range<-rev(c(c(0.8,1/2,1/4,1/10),rev(-c(0.8,1/2,1/4,1/10)))) ##proportional final change at end of interval
names(incidence.range)<-c(0.1,0.01,0.001,0.00001)
names(trend.range)<-trend.range

mean.incidence.range=c(0.1,0.01,0.001,0.00001)

ll.image.list<-vector(mode="list",length=length(incidence.range)*length(trend.range))
for(a in incidence.range[1]){
  for(b in trend.range[1]){
    list.ndx<-(match(a,incidence.range)-1)*length(trend.range)+match(b,trend.range)
    b.abs<-a*b/c

    
    list.ndx<-(match(a,incidence.range)-1)*length(trend.range)+match(b,trend.range)
    b.abs<-a*b/c
    
    plot.data<-na.omit(EndemicLinear(1000/(c*(a+b.abs*c/2)),
                                     incidence = a,
                                     trend = b.abs,start.time=0,end.time=c))
    
    a.range<-seq(1/2*a,1.5*a,length.out=res)
    b.range<-seq(max(-1/2*a/c,b.abs-abs(b.abs)/2),min(1/2*a/c,b.abs+abs(b.abs)/2),length.out=res)
    ll.grid<-outer(a.range,b.range,function(X,Y){
      sapply(1:length(X),function(x){
        tsi.linear.ll(pars =c(log(X[x]),log(1e-6+Y[x]+X[x]/c)),data = plot.data,censorLimit. = c,popsize =1000/(c*(a+b.abs*c/2)) )}
      )
      colnames(ll.grid)<-b.range
      rownames(ll.grid)<-a.range
      ll.image.list[[list.ndx]]<-ll.grid  
    })
  }
}

## Evaluating performance here




trendbias.grid<-outer(incidence.range,trend.range,function(X,Y){
  sapply(1:length(X),function(x){
    tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c,n.obs=3000)$trend.bias}
    )
  })

trendest.grid<-outer(incidence.range,trend.range,function(X,Y){
    sapply(1:length(X),function(x)tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c)$trend.est)
  })

incidence.est.grid<-outer(incidence.range,trend.range,function(X,Y){
  sapply(1:length(X),function(x)tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c)$incidence.est)
}
)

incidence.bias.grid<-outer(incidence.range,trend.range,function(X,Y){
  sapply(1:length(X),function(x)tsi.linear.est(X[x],X[x]*Y[x]/c,censorLimit=c)$incidence.bias)
}
)
est.df<-merge(melt(incidence.est.grid),melt(trendest.grid),by=c("Var1","Var2"))
names(est.df)<-c("true.inc","rel.trend","est.inc","est.trend")
est.df$true.trend<-est.df$true.inc*est.df$rel.trend/c
est.df$trend.relbias<-est.df$est.trend/est.df$true.trend
matplot(trend.range,t(trendbias.grid),type="l")

matplot(trend.range,t(trendest.grid),type="l")

melt(trendest.grid)
ggplot(melt(incidence.bias.grid),aes(y=value,x=Var2,group=Var1,col=factor(Var1)))+geom_line()
ggplot(melt(trendest.grid),aes(y=value/(Var2*Var1/c),x=Var2))+geom_line()+facet_wrap(~Var1,scales="free")
ggplot(melt(trendest.grid),aes(y=value,x=Var2,group=Var1,col=factor(Var1)))+geom_line()



###linear match

estTrend.df<- as.data.frame(apply(
  est.df,1,function(x)x[3]+x[4]*(0:10)))
colnames(estTrend.df)<-gsub("-","minus",paste(format(est.df$true.inc,scientific=FALSE,drop0trailing=TRUE),est.df$rel.trend,sep=".."))
estTrend.df<-data.frame(Time=0:10,estTrend.df)
estTrend.df<-melt(estTrend.df,id.vars="Time",variable.name="combo")
estTrend.df$true.inc<-as.numeric(gsub("X","",sapply(strsplit(as.character(estTrend.df$combo),"\\.\\."),"[",1)))
estTrend.df$rel.trend<-as.numeric(gsub("minus","-",sapply(strsplit(as.character(estTrend.df$combo),"\\.{2,}"),"[",2)))

ggplot(estTrend.df,aes(x=Time,y=value))+geom_line()+facet_grid(true.inc~rel.trend)+geom_abline(aes(intercept=true.inc,slope=true.inc*rel.trend/10),col="red")


##tsi match
estTrend.df<- as.data.frame(apply(
  est.df,1,function(x)lineardist(0:10,x[3],x[4])))
colnames(estTrend.df)<-gsub("-","minus",paste(format(est.df$true.inc,scientific=FALSE,drop0trailing=TRUE),est.df$rel.trend,sep=".."))
estTrend.df<-data.frame(Time=0:10,estTrend.df)
estTrend.df<-melt(estTrend.df,id.vars="Time",variable.name="combo")
estTrend.df$true.inc<-as.numeric(gsub("X","",sapply(strsplit(as.character(estTrend.df$combo),"\\.\\."),"[",1)))
estTrend.df$rel.trend<-as.numeric(gsub("minus","-",sapply(strsplit(as.character(estTrend.df$combo),"\\.{2,}"),"[",2)))
                                  
ggplot(estTrend.df,aes(x=Time,y=value))+geom_line()+facet_grid(true.inc~rel.trend)+
  geom_histogram(data=rlinear.tsi(1000,incidence=a,trend=b),aes(y=..count../tapply(..count..,..PANEL..,sum)[..PANEL..]),binwidth=1)
                                  
                                  
                                  

##Looking at plotting the likelihood function 

#This seems to indicate that the best estimate is biased and negative... 
test.data<-EndemicLinear(n.infections =10000,incidence = 0.01,trend=0.01/(30*2),start.time = 0,end.time = 30)
plot(seq(-0.01/30,0.01/30,length.out=100),
     sapply(seq(-0.01/30,0.01/30,length.out=100),
            function(x)tsi.linear.ll(c(log(0.01),log(x+2*0.01/30)),data=test.data,popsize=10000,censorLimit. = 30)
            )
     )

min.ndx<-which.min(sapply(seq(-0.01/30,0.01/30,length.out=100),
                 function(x)tsi.linear.ll(c(log(0.01),log(x+2*0.01/30)),data=test.data,popsize=10000,censorLimit. = 30)))
abline(v=0.01/(30*2),col="red")
abline(v=seq(-0.01/30,0.01/30,length.out=100)[min.ndx],col="blue")
seq(-0.01/30,0.01/30,length.out=100)[min.ndx]

hist(test.data,freq=FALSE,ylim=c(0,0.1))
for(beta in seq(-0.01/30,0.01/30,length.out=10)){lines(dlinear.tsi(InfTime =0:30,censorLimit = 30,incidence = 0.01,trend = beta),col="black",lwd=0.5)}
lines(dlinear.tsi(InfTime =0:30,censorLimit = 30,incidence = 0.01,trend = 0.005/30),col="red")
lines(dlinear.tsi(InfTime =0:30,censorLimit = 30,incidence = 0.01,trend = seq(-0.01/30,0.01/30,length.out=100)[min.ndx]),col="blue")
##Biased downards because there is no individuals infected at T==30...
