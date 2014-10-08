###Function by Gustaf Rydevik to generate simulated  population antibody
###and nucleic acid data,
### for some disease, including both infected and non-infected subjects.
### April-May 2011

LabdataGeneratorGeneric<-function(testkinetic,
                           timeFun,timeFun.args=NULL, 
                           errorFun,errorFun.args=NULL,...){

  #timeFun generates a set of infection times for a number of individuals

infection.times<-do.call("timeFun",c(timeFun.args,...))

##Generating expected mean test results
test.meanvalues<-testkinetic(infection.times)

test.obsvalues=do.call("errorFun",c(list(test.means.array=test.meanvalues),errorFun.args))
        
  return(list(infection.times=infection.times,
              test.meanvalues=test.meanvalues,
              test.obsvalues=test.obsvalues))
}
