
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
sim.dir<-"/Users/gustafrydevik/simulations/"
if(is.girion)sim.dir<-file.path(dropbox.path,"simulations")
##Project specific parameters
project.path<-file.path(dropbox.path,"PhD folder/Serology and NA")
if(is.girion){project.path<-dropbox.girion}
data.path<-file.path(project.path,"Data")
script.path<-file.path(project.path,"Code")
output.path<-file.path(project.path,"Output")

lapply(dir(file.path(dropbox.path,"GusLib"),full.names=T),source)
lapply(dir(file.path(script.path,"multitest library"),full.names=T),source)

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


# ##Scenario: 
# sample size 10-100, short,long, and full scenario, one test vs two tests. 
#

sim.scenarios<-list(pertussis.late=list(start.time=98,end.time=176,sample.size=c(10,20,40,80,89)),
                pertussis.late.full=list(start.time=98,end.time=252,sample.size=c(10,20,40,80,100,200,233)),
                btv.early=list(start.time=1,end.time=14,sample.size=c(10,20,24)),
                btv.full=list(start.time=1,end.time=49,sample.size=c(10,20,40,61))
  )
                
#sample.size.range<-c(seq(10,50,by=10),100,200,300)
sample.size.range<-c(20,50,100)

current.seed<-as.integer(as.numeric(format(Sys.time(),"%Y%m%d%H")))

    


for(rep in 0:19){
  for(simvars in sim.scenarios){
    for(Sample.size in sample.size.range){
      for(Onetest in c(1,0)){
        Pathogen<-if(simvars$start.time==1){"btv"}else{"pertussis"}
        rbatch(rfile=file.path(script.path,"SerologyDataSim_lognorm_fixedvariance.R"),
               seed=current.seed,
               samplesize=Sample.size,
               Actual.=1,
               onetest=Onetest,
               nreps=5,
               rep.prefix=rep,
               Start.time=simvars$start.time,
               Set.endtime=simvars$end.time,
               burn.in=8000,
               adapt.iter=2000,
               converge=F,
               converge.criteria=1.15,
               converge.time=5,
               pathogen=Pathogen)
      }
    }
  }
}
rbatch.local.run()
