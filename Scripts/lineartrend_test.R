library(IDPmisc)
###from 
#http://theoreticalecology.wordpress.com/2011/12/09/mcmc-chain-analysis-and-convergence-diagnostics-with-coda-in-r/

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="blue4", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = "spearman"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r)
}

betterPairs <- function(YourData){
  return(pairs(YourData, lower.panel=function(...) {par(new=TRUE);ipanel.smooth(...)}, diag.panel=panel.hist, upper.panel=panel.cor))
}

# Example

x = rnorm(10000)
betterPairs(data.frame(A = x, B = 0.6 * x + 0.3 * rnorm(10000), C = rnorm(10000)))

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
sim.dir<-"/Users/gustafrydevik/simulations/"
if(is.girion)sim.dir<-file.path(dropbox.path,"simulations/")
##Project specific parameters
project.path<-file.path(dropbox.path,"PhD folder/Chapter3")
if(is.girion){project.path<-dropbox.girion}
data.path<-file.path(project.path,"Data")
script.path<-file.path(project.path,"Scripts")
output.path<-file.path(project.path,"Output")

lapply(dir(file.path(dropbox.path,"GusLib"),full.names=T),source)


autolib(rjags)
autolib(batch)
autolib(reshape2)

My.device<-Gen.device("png",res=400,width=12,height=3,units="in")



kinetic.fun<-LotkaVolterra.Fun()
end.time<-100
incidence=1/100
pop<-1:1e5
InfTime<-rep(1e5,length(pop))

iteration.data<-
  LabdataGeneratorGeneric(
    testkinetic=kinetic.fun,
    timeFun=EndemicLinear,
    timeFun.args=list(n.infection=1000,
                      start.time=1,
                      end.time=30,
                      incidence=10/100,
                      trend=0
                      ),
    errorFun=errorFun.lognorm,
    errorFun.args=list(standard.deviation=log(c(1.2,1.1)))#log(measurement.sd))
  )
bugsdata<-list(N=nrow(iteration.data$test.obsvalues),
               Test.data=iteration.data$test.obsvalues,
               time.lookup=seq(0.5,99.5,by=0.5),
               test.lookup=kinetic.fun(c(seq(0.5,99.5,by=0.5))),
               is.naive=is.na(iteration.data$test.obsvalues[,1]),
               ntest=n.tests,censorLimit=30)
pars.inits<-vector(length=n.chains.,mode="list")
for(C in 1:n.chains.){
  pars.inits[[C]]<-list(
    sd.exp.tmp=c(1.2,1.1),#exp(runif(n.tests,0.02,2)),
    InfTime=ifelse(is.na(iteration.data$test.obsvalues[,1]),
                   35,runif(nrow(iteration.data$test.obsvalues),0,30)),
    incidence=runif(1,0,1),
    trend=runif(1,-0.01/35,0.01/35)
  )
}


tmp<-jags.model(file=file.path(script.path,"bugs/hindcast_linear_centered.txt"),
           data=bugsdata,inits=pars.inits,
           n.adapt=300,
           n.chains=5)


###Seems to take ~5-10k iterations w. the hindcast_linear_centered.txt file and normal prior for incidence. 
tmp.samples<-coda.samples(tmp,variable.names=c("incidence","trend","InfTime","sd","tau"),n.iter=1000,thin=5)

lineardist<-function(time,alpha,beta,chain){
  alpha<-alpha[chain]
  beta<-beta[chain]
  (alpha+beta*time)*exp(-(alpha+beta*time/2)*time)}



##diagnostics

par(mfrow=c(1,2))
for(i in 1:5){
  betterPairs(as.data.frame(tmp.samples[[i]][,c("incidence","trend")]))
}
par(mfrow=c(2,3))
for(i in 1:5){
  hist(iteration.data$infection.times,freq=FALSE)
  #hist(tmp.samples[[i]][,grep("InfTime",colnames(tmp.samples[[1]]))],freq=FALSE)
  lapply(lapply(tmp.samples[,c("incidence","trend")],colMeans),function(x)
    lines(seq(1,30,by=1),lineardist(seq(1,30,by=1),x[1],-x[2],1),lwd=2,col="red"))
  
}
par(mfrow=c(1,1))
plot(tmp.samples[,c("incidence","trend","sd[1]","sd[2]")])

par(mfrow=c(2,3))
for(i in 1:5){
  ###observed test measuremnts vs estimated means
  plot(iteration.data$test.obsvalues[,1],
       kinetic.fun(colMeans(tmp.samples[[i]][,grep("InfTime",colnames(tmp.samples[[5]]))]))[,1])
}
##infection times distribution


par(mfrow=c(1,2))
for(i in 1:5){
##True vs estimated infeciton times, and true infection time vs error in infection time estimate
plot(iteration.data$infection.times,colMeans(tmp.samples[[i]][,grep("InfTime",colnames(tmp.samples[[1]]))]))
plot(iteration.data$infection.times,iteration.data$infection.times-colMeans(tmp.samples[[i]][,grep("InfTime",colnames(tmp.samples[[5]]))]))
}


tmp<-jags.model(file.path(script.path,"bugs/hindcast_linear_test.txt"),
                data=list(N=length(increase.true.inftime[,1]),censorLimit=30,
                          InfTime=increase.true.inftime[,1]),
                inits=list(trend.tmp=1/100/30,incidence=1/100))
adapt(tmp,5000)
update(tmp,10000)
tmp.samples<-coda.samples(tmp,variable.names=c("incidence","trend","L"),50000)
hist(increase.true.inftime[,1],breaks=seq(0,31,by=1),xlim=c(0,100),freq=F)
est.dist<-lineardist(1:100,mean(tmp.samples[[1]][,"incidence"]),mean(tmp.samples[[1]][,"trend"]),chain=1)
lines(1:100,lineardist(1:100,mean(tmp.samples[[1]][,"incidence"]),mean(tmp.samples[[1]][,"trend"]),chain=1)/sum(est.dist[1:30]))

lines(1:100,lineardist(1:100,1/100,1/100/50,chain=1)*1000/329,col="red")


linear.ll<-function(incidence,trend,InfTime){
  censorLimit=30
  dunif(incidence,log=T)+
    dnorm(trend,0,1/sqrt(0.0001),log=T)+
  log((incidence-trend*InfTime)*exp(-(incidence-trend*(InfTime)/2)*InfTime)/(1-exp(-(incidence+trend*censorLimit/2)*censorLimit)))
}

linear.ll(1/100,-1/100/50,InfTime=increase.true.inftime[,1])
linear.ll(mean(tmp.samples[[1]][,1]),-mean(tmp.samples[[1]][,2]),InfTime=increase.true.inftime[,1])

InfTime<-increase.true.inftime[,1]
incidence<-1/100
trend<--1/100/50
(incidence-trend*InfTime)*exp(-(incidence-trend*(InfTime)/2)*InfTime)
InfTime<-sort(InfTime)
InfTime
linear.ll(1/100,-1/100/50,InfTime=increase.true.inftime[,1])
(incidence-trend*InfTime)*exp(-(incidence-trend*(InfTime)/2)*InfTime)
incidence<-mean(tmp.samples[[1]][,1])
trend<-mean(tmp.samples[[1]][,2])
(incidence-trend*InfTime)*exp(-(incidence-trend*(InfTime)/2)*InfTime)


