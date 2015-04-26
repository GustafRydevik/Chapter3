
####Use 
#bash > nohup R --vanilla --slave --args seed 1000 burn.in 25000 adapt.iter 500 Set.endtime 63 nreps 200 < SerologyDataSim_Pertussis.R > nohup2.out 2>&1&
##to run in girion

###parameters

####Generating various simulated testing datasets, using functions from eventGenerators.R and Labdatagenerator.

##Parameters for using the batch package for multiple parallell runs


###The following scenarios I want to run:


# 1)  The simple increasing, decreasing and constant example
# 
# Repeats (x10-x100 each): 
#   - All have same level of noise
# -
#   
#   1) Effect of sample size by slope and diseasetype
# 100,250,500,1000, 2000 samples
# +/- 5%,2.5%,0% slope 
# 
# 5x3*5=125
# 10 reps
# 
library(data.table)
sample.size<-c(100,250,500,1000,2000)
slopes<-c(-5,-2.5,0,2.5,5)/10
diseaseType=c("diseaseType1","diseaseType2","diseaseType3")
scenarios.pars.slopeDiseaseType.df<-data.frame(expand.grid(sample.size=sample.size,slope=slopes,diseaseType=diseaseType),
                              ntest="Both",
                              duration=10,
                              batch.name="SlopeDiseaseType",
                              reps=10,
                              intercept=0.1)


# 2)  effect of disease type and number of tests by sample size, with fixed slope
# 100,250,500,1000, 2000 samples
# 
# 3*3*5=45
# 10 reps
# 

sample.size<-c(100,250,500,1000,2000)
ntest<-c("Ab","NA","Both")
diseaseType=c("diseaseType1","diseaseType2","diseaseType3")
scenarios.pars.ntestDiseaseType.df<-data.frame(expand.grid(sample.size=sample.size,ntest=ntest,diseaseType=diseaseType),
                                               duration=10,
                                               batch.name="ntestDiseaseType",
                                               reps=10,
                                               slope=-0.5,
                                               intercept=0.1)



# 3) effect of length of duration and type of disease
# 10,15,20 units 
# 100,250,500,1000, 2000 samples
# 3*3*5=45
# 
# 
sample.size<-c(100,250,500,1000,2000)
ntest<-c("Ab","NA","Both")
diseaseType=c("diseaseType1","diseaseType2","diseaseType3")
duration<-c(10,15,20)
scenarios.pars.durationDiseaseType.df<-data.frame(expand.grid(sample.size=sample.size,duration=duration,diseaseType=diseaseType),
                                               ntest="Both",
                                               batch.name="durationDiseaseType",
                                               reps=10,
                                               slope=-0.5,
                                               intercept=0.1)



##4) "realistic scenarios" : Squirrel pox, Scrapie and Chlamydia 

##Squirrel pox

sample.size<-c(100,250)
slope<-c(-0.5,-0.25,-0.1,0,0.1,0.25,0.5)
intercept=0.3
diseaseType=c("diseaseTypeSquirrelpox")##Code this disease! 
duration<-c(3)

scenarios.pars.squirrelpox.df<-data.frame(expand.grid(sample.size=sample.size,slope=slope,duration=duration,diseaseType=diseaseType),
                                                  ntest="Both",
                                                  batch.name="Squirrelpox",
                                                  reps=1,
                                                  intercept=0.3)

##Chlamydia

sample.size<-c(100,250,500,1000,2000,5000,10000)
slope<-c(0.5)
intercept=0.02
diseaseType=c("diseaseTypeChlamydia")##Code this disease! 
duration<-c(10)

scenarios.pars.chlamydia.df<-data.frame(expand.grid(sample.size=sample.size,slope=slope,duration=duration,diseaseType=diseaseType),
                                          ntest="Both",
                                          batch.name="Chlamydia",
                                          reps=1,
                                          intercept=0.02)
##Scrapie

sample.size<-c(100,250,500,1000,2000,5000,10000)
slope<-c(-0.75)
intercept=0.05/100
diseaseType=c("diseaseTypeScrapie")##Code this disease! 
duration<-c(6)

scenarios.pars.scrapie.df<-data.frame(expand.grid(sample.size=sample.size,slope=slope,duration=duration,diseaseType=diseaseType),
                                        ntest="Both",
                                        batch.name="Scrapie",
                                        reps=1,
                                        intercept=0.05/100)





scenario.pars<-do.call("rbind",list(a=scenarios.pars.slopeDiseaseType.df,b=scenarios.pars.ntestDiseaseType.df,c=scenarios.pars.durationDiseaseType.df,
                                    d=scenarios.pars.squirrelpox.df,
                                    e=scenarios.pars.chlamydia.df,
                                    f=scenarios.pars.scrapie.df))
scenario.pars<-data.table(scenario.pars)


# (all ss together is ~ 4000 individuals, or 1/2 of a full 10k run - thus, this would equivalent to running 10*(15+9+9)/2=330/2=165 10k runs.  Comparatively, we ran 322 complete runs last time, in addition to all the smaller scenarios.
#  
#  In addition, we would like to have  a 2-3 of ”realistic” scenarios, with fixed incidence, fixed slope, fixed test, but varying sample size. However, this would only be a small addition in computation time
#  
#  
#  The structure would thus be 
#  
#  1) simple results 
# 2) example from specific scenarios inspired by reality
# 3) Effect of type of test kinetic and number of tests used 
# 4) How far back  can we hindcast depending on test kinetic?
# 5) Effect of size of trend

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


base.pars<-list(
  ### MCMC parameters
  burn.in=10,
  adapt.iter=10,
  n.chains.=1,
  mcmc.ss=10,
  thin=5,
  
  ###Other run parameters
  converge=FALSE,
  converge.criteria=1.15,
  converge.time=0.2
)
scenario.pars$sample.size<-
  scenario.pars$sample.size/10

scenario.pars$sample.size<-
  scenario.pars$sample.size*10
base.pars<-list(
  ### MCMC parameters
  burn.in=100,
  adapt.iter=100,
  n.chains.=4,
  mcmc.ss=1000,
  thin=5,
  
  ###Other run parameters
  converge=FALSE,
  converge.criteria=1.15,
  converge.time=0.2
)

seed<-1000
seed.iter<-seed
Start.time=0

Number.cores<-4

if(is.amz){Number.cores<-as.numeric(system("nproc",intern=TRUE))}

library(data.table)
scenario.pars<-as.data.table(scenario.pars)
seed<-1000
for(row in  1:nrow(scenario.pars)){
               scenario.args<-list(
               Start.time=0,
               End.time=scenario.pars[row,duration],
               Incidence=scenario.pars[row,intercept],
               change.percent=scenario.pars[row,slope],
               Trend=scenario.pars[row,intercept*slope/(duration)],
               diseaseType=as.character(scenario.pars[row,diseaseType]),
               nreps=scenario.pars[row,reps],
               Epi.scenario="EndemicLinear")
               all.args<-c(list(rfile=shQuote(file.path(script.path,"hindcasting_generic_awsmarch.R"))),
                           base.pars,
                           samplesize=scenario.pars[row,sample.size],
                           scenario.args,
                           list(
                             nreps=scenario.pars[row,reps],
                             ncores=4),
                           list(test1.sd=1.25,
                                test2.sd=1.25,
                                n.tests=scenario.pars[row,ntest]),
                           seed=seed.iter
               )
          do.call("rbatch",all.args)
                  
          
          seed.iter<-seed.iter+1
        }

rbatch.local.run(ncores=Number.cores)
