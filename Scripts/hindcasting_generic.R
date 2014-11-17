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
sim.dir<-"/Users/gustafrydevik/simulations/Chapter3/"
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

source(file.path(script.path,"hindcasting helper functions/DensityEstimator.R"))
#### Reading in  functions specific for this project 
lapply((grep("R$",dir(file.path(script.path,"hindcasting\ helper\ functions"),full.names=T),value=T)),source)


########################
###Default parameters ##
########################

### MCMC parameters 
seed=1000
adapt.iter=1000
n.chains.=5
mcmc.ss=1000
thin=1

###Other run parameters
converge<-F
converge.criteria<-1.15
converge.time=0.2
rep.prefix=NULL
nreps=0

###Epidemic trend pars
Epi.scenario="EndemicConstantFixedSD" ##EndemicLinear EpidemicExp EpidemicLognorm
Epi.model="EndemicLinear"
Start.time=1
End.time=30
Incidence=1/100
samplesize= 1000
##Par for linear and exponential trend
Trend=0.01/50
Peak.time=15
Epi.sd=10



##Pars releveant for trends modelled on true outbreaks
Actual.=T
pathogen="btv" #or "pertussis"

##Test value generator parameters 
diseaseType="diseaseType1" ##diseaseType2 diseaseType3
n.tests=2
test1.sd=1.27
test2.sd=1.57




###Alternatively, reading in all parameters from a batch call (how many pars can you pass?)
parseCommandArgs()
set.seed(seed)



##LV.Fun accepts a list of parameters specific to the disease type. get() gets this list.
kinetic.fun=LotkaVolterra.Fun(disease=get(diseaseType))
measurement.sd=c(test1.sd=1.27,test2.sd=1.57)


####################################################################
######## setting parameters for the different scenarios ############
####################################################################



#Constant trend args

if(Epi.scenario%in%c("EndemicConstant","EndemicIncrease","EndemicDecrease")){
##Linear trend args
scenario.args=list(n.infection=1000,
                   start.time=Start.time,
                   end.time=End.time,
                   incidence=Incidence,
                   trend=Trend)
modelscript.name="bugs/hindcast_linear.txt"
sample.vars=c("incidence","trend","sd","InfTime")
Epi.model="EndemicLinear"
}else{
  if(Epi.scenario%in%c("EndemicConstantFixedSD","EndemicIncreaseFixedSD","EndemicDecreaseFixedSD")){
    ##Linear trend args
    scenario.args=list(n.infection=1000,
                       start.time=Start.time,
                       end.time=End.time,
                       incidence=Incidence,
                       trend=Trend)
    modelscript.name="bugs/hindcast_linear_fixedSD.txt"
    sample.vars=c("incidence","trend","InfTime")
    Epi.model="EndemicLinear"
  }else{
##Exponential trend args
  if(Epi.scenario=="EpidemicExp"){
scenario.args=list(n.infection=samplesize,
                   start.time=Start.time,
                   end.time=End.time,
                   trend=Trend)
modelscript.name="bugs/hindcast_exponential.txt"
sample.vars=c("trend")
Epi.model="EpidemicExp"

}else{
##lognorm trend args 
scenario.args=list(n.infection=samplesize,
                   start.time=Start.time,
                   end.time=End.time,
                   peak.time=Peak.time,
                   sd=Epi.sd)
                 
Epi.model="EpidemicLognorm"

modelscript.name="bugs/hindcast_lognorm.txt"
sample.vars=c("peak.time","duration")
}}}
##############################

#### Step.data (and likelihood!)


scenarios.results<-vector(mode="list",length=0)


for(i in 0:nreps){
 
  t0<-proc.time()
        iteration.data<-
          LabdataGeneratorGeneric(
            testkinetic=kinetic.fun,
            timeFun= get(Epi.model),
            timeFun.args=scenario.args,
            errorFun=errorFun.lognorm,
            errorFun.args=list(standard.deviation=log(measurement.sd))
          )
  n.sampled<-length(iteration.data$infection.times)
  iteration.data<-lapply(iteration.data,function(x){if(!is.null(ncol(x))){x[is.finite(iteration.data$infection.times),]}else{x[is.finite(iteration.data$infection.times)]}})
  iteration.data$n.sampled<-n.sampled
  iteration.data$PropAboveCensor<-(n.sampled-length(iteration.data$infection.times))/n.sampled
  ###########################################################
  ##### Generating inits and bugs pars below here ###########
  ###########################################################
  
 
   if(Epi.scenario%in%c("EndemicConstant","EndemicIncrease","EndemicDecrease")){
      ##Linear trend args
      bugs.args=list(censorLimit=End.time-Start.time+1,
                     is.naive=is.na(iteration.data$test.obsvalues[,1])
      )
      inits.list<-vector(mode="list",length=n.chains.)
      for(C in 1:n.chains.){
      inits.list[[C]]=list(
        InfTime=ifelse(is.na(iteration.data$test.obsvalues[,1]),
                       35,runif(nrow(iteration.data$test.obsvalues),0,30)),
        incidence=runif(1,0.01,0.3),
        trend=runif(1,-0.01/35,0.01/35))
      }
    }else{
      if(Epi.scenario%in%c("EndemicConstantFixedSD","EndemicIncreaseFixedSD","EndemicDecreaseFixedSD")){
        ##Linear trend args
        bugs.args=list(censorLimit=End.time-Start.time+1,
                       is.naive=is.na(iteration.data$test.obsvalues[,1]),
                       sd=c(test1.sd,test2.sd)
        )
        inits.list<-vector(mode="list",length=n.chains.)
        for(C in 1:n.chains.){
          inits.list[[C]]=list(
            InfTime=runif(nrow(iteration.data$test.obsvalues),0,End.time-Start.time+1),
            incidence=runif(1,0.01,0.3),
            trend=runif(1,-0.01/35,0.01/35))
        }
      }else{
      ##Exponential trend args
      if(Epi.scenario=="EpidemicExp"){
        bugs.args=list()
        inits.list=list()

        
      }else{
        ##lognorm trend args 
        bugs.args=list()
        inits.list=list(peak.time=rgamma(1,4,scale=End.time/4),
                        duration.tmp=abs(rt(1,5))
        )
        
      }}}
  ##########################################
  ##########################################
          bugsdata<-c(list(N=nrow(iteration.data$test.obsvalues),
                 Test.data=iteration.data$test.obsvalues,
                 time.lookup=c(seq(0.5,99.5,by=0.5)),
                 test.lookup=kinetic.fun(c(seq(0.5,99.5,by=0.5))),
                 ntest=n.tests),
                 bugs.args
                )
        pars.inits<-vector(length=n.chains.,mode="list")
        for(C in 1:n.chains.){
          pars.inits[[C]]<-c(inits.list[[C]],
                             list(logsd.tmp=rep(runif(1,log(1.02),log(4)),n.tests))
          )
        }
 
        
          multitest.bugs<-jags.model(file=file.path(script.path,modelscript.name),
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
          iteration.samples<-coda.samples(multitest.bugs,variable.names=sample.vars,mcmc.ss,thin=thin)
        }
       gelman.current<-mean(gelman.diag(iteration.samples[,grep(paste(sample.vars[-which(sample.vars=="InfTime")],collapse="|"),colnames(iteration.samples[[1]])),drop=FALSE])$psrf[,1])
            print("current convergence: ")
        print(gelman.current)
        print("\n")
        if(converge==F){gelman.current<-1}
        elapsed.time<-(proc.time()-start.time)[3]/3600 
        print(elapsed.time)
        }
          scenarios.results$true.values<-iteration.data
          scenarios.results$est<-iteration.samples
          
        
        scenarios.results$pars.inits<-pars.inits
        if((inherits(possibleError1, "error")|inherits(possibleError2, "error"))){
          scenarios.results$est<-FALSE
          scenarios.results$error<-c(possibleError1,possibleError2)
        }  
        
      
      print(i) 
      save(scenarios.results,file=paste(
        sim.dir,diseaseType,
        "_",Epi.scenario,
          "_seed_",seed,
          "_t",n.tests,
          paste("_samplesize",samplesize,sep=""),
          "_",rep.prefix,i,
        "_gr",round(gelman.current,2),
          ".RData",sep=""))
      }
