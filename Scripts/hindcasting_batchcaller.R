
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



base.pars<-list(
### MCMC parameters
burn.in=1000,
adapt.iter=500,
n.chains.=5,
samplesize= 10000,
mcmc.ss=10000,

###Other run parameters
converge=FALSE,
converge.criteria=1.15,
converge.time=0.2,


###Epidemic trend pars
Start.time=1,
End.time=30,
Incidence=1/100
)

constant.pars<-list(
Epi.scenario="EndemicConstant" ##EndemicLinear EpidemicExp EpidemicLognorm
Trend=0
)
##Par for linear increase
linear.increase.pars<-list(
  Epi.scenario="EndemicIncrease", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=-1/100/50
)

##Par for linear decrease
linear.decrease.pars<-list(
  Epi.scenario="EndemicDecrease", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=1/100/50
)


constant.pars.fixedSD<-list(
  Epi.scenario="EndemicConstant", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=0
)
##Par for linear increase
linear.increase.pars.fixedSD<-list(
  Epi.scenario="EndemicIncrease", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=-1/100/50
)

##Par for linear decrease
linear.decrease.pars.fixedSD<-list(
  Epi.scenario="EndemicDecrease", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=1/100/50
)

##Exponential trend
exponential.pars<-list(
  Epi.scenario="EpidemicExp", ##EndemicLinear EpidemicExp EpidemicLognorm
  Trend=2
)
##lognormal trend
lognormal.pars<-list(
  Epi.scenario="EpidemicLognorm", ##EndemicLinear EpidemicExp EpidemicLognorm
  Peak.time=15,
  Epi.sd=10
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
test1.sd=1.27
test2.sd=1.57

sample.size.range<-c(10000,20000,50000)[1]
ntest.range=2
    

reps.per.call=1
ncalls.per.combination=1
seed.iter<-seed
for(rep in (ncalls.per.combination-1)){
  for(Sample.size in sample.size.range){
    for(trendpars in list(constant.pars.fixedSD,linear.increase.pars.fixedSD,linear.decrease.pars.fixedSD)){
      for(diseasepars in list(d1.pars,d2.pars,d3.pars)[1]){
        for(ntests in ntest.range[1]){
          do.call("rbatch",
                  c(list(rfile=shQuote(file.path(script.path,"hindcasting_generic.R"))),
                    base.pars,
                    trendpars,
                    diseasepars,
                    list(rep.prefix=rep,
                         nreps=(reps.per.call-1),
                         ncores=1),
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
}
rbatch.local.run(ncores=4)
