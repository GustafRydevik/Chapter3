##Accepts a vector and/or an array, as well as a vector of standard deviations that is recycled if it is shorter than the number of columns

errorFun.lognorm<-function(test.means.array,standard.deviation=log(1.1)){
  if(is.null(ncol(test.means.array))){
    test.means.array<-array(test.means.array,dim=c(length(test.means.array),1))
  }
  col=1
  obs.value.array<-apply(test.means.array,2,function(x){
    obs.error=rlnorm(length(x),meanlog=0,sdlog=standard.deviation[col])
    obs.value=obs.error*x
    col<<-ifelse((col+1)>length(standard.deviation),1,(col+1))
    return(obs.value)
  })
  return(obs.value.array)
}
