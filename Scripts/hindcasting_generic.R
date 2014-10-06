####Use 
#bash > nohup R --vanilla --slave --args seed 1000 burn.in 25000 adapt.iter 500 Set.endtime 63 nreps 200 < SerologyDataSim_btv.R > nohup2.out 2>&1&
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
sim.dir<-"/Users/gustafrydevik/simulations/"
if(is.girion)sim.dir<-file.path(dropbox.path,"simulations/")
##Project specific parameters
project.path<-file.path(dropbox.path,"PhD folder/Chapter3")
if(is.girion){project.path<-dropbox.girion}
data.path<-file.path(project.path,"Data")
script.path<-file.path(project.path,"Scripts")
output.path<-file.path(project.path,"Output")

lapply(dir(file.path(dropbox.path,"GusLib"),full.names=T),source)

##Note to self: pare this down substantially!
# autolib(mclust)
# autolib(lattice)
# autolib(KernSmooth)
# autolib(RColorBrewer)
# autolib(deldir)
# autolib(spatstat)
# autolib(maptools)
# autolib(sp)
autolib(rjags)
autolib(batch)
# autolib(ggplot2)
# autolib(gridExtra)
autolib(reshape2)

My.device<-Gen.device("png",res=400,width=12,height=3,units="in")

#### Reading in  functions specific for this project 
lapply(grep("R$",dir(file.path(script.path,"hindcasting helper functions/"),full.names=T),value=T),source)


seed<-1000
nreps=1
ntests=2
measurement.sd<-c(t1=1.27,t2=1.57)

burn.in<-1000 ##10000 seems enough for endtime=63
adapt.iter<-1000
pathogen="btv" #or "pertussis"

End.time=100#7*7
Start.time=1
n.chains.=5
samplesize= 1000
Actual.=T
n.tests=2
converge<-T
converge.criteria<-1.15
converge.time=0.2
rep.prefix=NULL



kinetic.fun<-LotkaVolterra.Fun(disease=diseaseType1)
epidemic.trend.fun<-pertussis.wi.Fun

parseCommandArgs()
set.seed(seed)

#### Step.data (and likelihood!)

scenarios.data<-vector(mode="list",
                           length=1)

scenarios.results<-vector(mode="list",length=1)


for(i in 0:nreps){
 
  t0<-proc.time()
        iteration.data<-
          LabdataGeneratorGeneric(
            testkinetic=kinetic.fun,
            timeFun=EndemicDisease,
            timeFun.args=list(n.infection=samplesize,
                              start.time=Start.time,
                              end.time=Set.endtime,
                              Actual=Actual.,
                              incidence=1/100),
            errorFun=errorFun.lognorm,
            errorFun.args=list(standard.deviation=log(c(1.01,1.01)))#log(measurement.sd))
          )
          bugsdata<-list(N=nrow(iteration.data$test.obsvalues),
                 Test.data=iteration.data$test.obsvalues,
                 time.lookup=c(seq(0.5,99.5,by=0.5),seq(100,198,by=2),Inf),
                 test.lookup=kinetic.fun(c(seq(0.5,99.5,by=0.5),Inf)),
                 prior.start=2^round(log2(Set.endtime-Start.time)),
                 ntest=n.tests)
        pars.inits<-vector(length=n.chains.,mode="list")
        for(C in 1:n.chains.){
          pars.inits[[C]]<-list(
            peak.time=rgamma(1,4,scale=Set.endtime/4),
            duration.tmp=abs(rt(1,5)),
            logsd.tmp=rep(runif(1,log(1.02),log(4)),n.tests))
        }
 
        
          multitest.bugs<-jags.model(file=file.path(script.path,"bugs/hindcast_constant.txt"),
                                              data=bugsdata,
                                              inits=pars.inits,
                                              n.adapt=0,
                                              n.chains=n.chains.)
         
        possibleError1<-tryCatch(adapt(multitest.bugs,adapt.iter),
                                 error=function(e) e )
        
        gelman.current<-Inf
        start.time<-proc.time()
        elapsed.time<-0
        while(((gelman.current>converge.criteria)&(elapsed.time<converge.time))){
        
        possibleError2<-tryCatch(update(multitest.bugs,burn.in),
                                 error=function(e) e)
        if(!(inherits(possibleError1, "error")|inherits(possibleError2, "error"))){
          iteration.samples<-coda.samples(multitest.bugs,variable.names=mod.varnames,1000)
        }
        }
       gelman.current<-gelman.diag(iteration.samples[,c("EpiStart","lambda"),drop=FALSE])$mpsrf
            print("current convergence: ")
        print(gelman.current)
        print("\n")
        if(converge==F){gelman.current<-1}
        elapsed.time<-(proc.time()-start.time)[3]/3600 
        print(elapsed.time)
        }
          scenarios.results[[as.character(SD)]]$true.values<-iteration.data
          scenarios.results[[as.character(SD)]]$est<-iteration.samples
          
        
        scenarios.results[[as.character(SD)]]$pars.inits<-pars.inits
        if((inherits(possibleError1, "error")|inherits(possibleError2, "error"))){
          scenarios.results[[as.character(SD)]]$est<-FALSE
          scenarios.results[[as.character(SD)]]$error<-c(possibleError1,possibleError2)
        }  
        
      
      print(i) 
      save(scenarios.results,file=paste(
        sim.dir,pathogen,
          "_seed_",seed,
          ifelse(onetest,"_one","_two"),
          "test_starttime_",Start.time,
          "_endtime_",Set.endtime,
          paste("_samplesize",samplesize,sep=""),
          "_",rep.prefix,i,
          ".RData",sep=""))
      }
