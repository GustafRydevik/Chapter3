
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
is.amz<-grepl("^ip",system("hostname",intern=TRUE))
dropbox.amz<-"/home/gustaf/ch3"
dropbox.sac<-"C:/Documents and Settings/GRydevik/My Documents/Dropbox"
dropbox.bioss<-"D:\\Dropbox"
dropbox.osx<-"/Users/gustafrydevik/Dropbox"
dropbox.girion<-"/home/gustaf/SeroNA_temp"
dropbox.path<-c(dropbox.osx,dropbox.sac,dropbox.bioss)[got.mdrive+is.win+1]
if(is.girion){dropbox.path<-dropbox.girion}
if(is.amz){dropbox.path<-dropbox.amz}
#sim.dir<-"D:/simulations/"
#sim.dir<-paste(dropbox.path,"/simulations",sep="")
sim.dir<-"/Users/gustafrydevik/simulations/Chapter3/"
if(is.girion)sim.dir<-file.path(dropbox.path,"simulations/")
if(is.amz)sim.dir<-file.path(dropbox.path,"simulations/")

##Project specific parameters
project.path<-file.path(dropbox.path,"PhD folder/Chapter3")
if(is.girion){project.path<-dropbox.girion}
if(is.amz){project.path<-file.path(dropbox.amz,"Chapter3")}
data.path<-file.path(project.path,"Data")
script.path<-file.path(project.path,"Scripts")
output.path<-file.path(project.path,"Output")

#lapply(dir(file.path(dropbox.path,"GusLib"),full.names=T),source)
autolib<-function(Package,mirror="http://stat.ethz.ch/CRAN",...){
  if(!suppressWarnings(require(as.character(substitute(Package)),
                               character.only=TRUE,quietly=TRUE,...))){
    Package<-as.character(substitute(Package))
    install.packages(Package,repos=mirror,...)
    library(Package,character.only=TRUE)	
  }		
}


lapply(dir(file.path(script.path,"hindcasting helper functions"),full.names=T),source)

##Note to self: pare this down substantially!
autolib(rjags)
autolib(batch)
autolib(reshape2)

#My.device<-Gen.device("png",res=400,width=12,height=3,units="in")

#### Reading in  functions specific for this project 


constant.pars<-list(
Epi.scenario="EndemicConstant", ##EndemicLinear EpidemicExp EpidemicLognorm
Trend=0
)
##Par for linear increase
linear.increase.pars<-list(
  Epi.scenario="EndemicIncrease", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=-5/100/25
)

##Par for linear decrease
linear.decrease.pars<-list(
  Epi.scenario="EndemicDecrease", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=5/100/25 #Trend as measured by time Now -> Past
)


constant.pars.fixedSD<-list(
  Epi.scenario="EndemicConstantFixedSD", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=0
)
##Par for linear increase
linear.increase.pars.fixedSD<-list(
  Epi.scenario="EndemicIncreaseFixedSD", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend= -5/100/30 #Trend as measured by Time Now -> Past 
)

##Par for linear decrease
linear.decrease.pars.fixedSD<-list(
  Epi.scenario="EndemicDecreaseFixedSD", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=5/100/30
)




##Test value generator parameters

d1.pars<-list(
diseaseType="diseaseType1" ##diseaseType2 diseaseType3
)

d2.pars<-list(
  diseaseType="diseaseType3" ##diseaseType2 diseaseType3
)

d3.pars<-list(
  diseaseType="diseaseType3" ##diseaseType2 diseaseType3
)

##test kinetics pars
 ##Fix this so we can refer to either NA or AB




base.pars<-list(
  ### MCMC parameters
  burn.in=100,
  adapt.iter=1000,
  n.chains.=5,
  mcmc.ss=1000,
  thin=5,
  
  ###Other run parameters
  converge=FALSE,
  converge.criteria=1.15,
  converge.time=0.2
)

ntest.range=c("Ab","NA","Both")
sample.size.range<-c(100,1000,10000)
change.percent.range<-c(-50,-25,0,25,50)/100 ##trend percentage change
intercept.par.range<-c(0.1,0.01,0.001)

End.time.range<-c(10,15,20)
obs.sd.range<-c(1.05,1.25,1.5)



reps.per.call=1
ncalls.per.combination=1
seed<-1000
seed.iter<-seed
Start.time=0

Number.cores<-4

if(is.amz){Number.cores<-as.numeric(system("nproc",intern=TRUE))}

#### 1215 combinations are too many!!!!!
for(rep in (ncalls.per.combination-1)){
  for(Sample.size in sample.size.range){
    for(change.percent in change.percent.range){
      for(intercept.par in intercept.par.range){
        for(diseasepars in list(d1.pars,d2.pars,d3.pars)){
          for(End.time in End.time.range){
             for(ntests in ntest.range){
               scenario.args<-list(start.time=Start.time,
                                   Start.time=Start.time,
               End.time=End.time,
               Incidence=intercept.par,
               Trend=intercept.par*change.percent/(End.time-Start.time),
               Epi.scenario="EndemicLinear")
               
          do.call("rbatch",
                  c(list(rfile=shQuote(file.path(script.path,"hindcasting_generic.R"))),
                    base.pars,
                    diseasepars,
                    samplesize=Sample.size,
                    scenario.args,
                    change.percent=change.percent,
                    list(rep.prefix=rep,
                         nreps=(reps.per.call-1),
                         ncores=8),
                    list(test1.sd=test1.sd,
                         test2.sd=test2.sd,
                         n.tests=ntests),
                    seed=seed.iter
                  )
          )
          seed.iter<-seed.iter+1
        }
      }
    }
  }
}}}
rbatch.local.run(ncores=Number.cores)
