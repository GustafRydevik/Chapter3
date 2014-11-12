
###Lotka-Volterra example kinetics
require(deSolve)
time.scale<-7
 Pars1<-c(growth.na=10/time.scale,killed.na=0.5/time.scale,feed.ab=2/time.scale,dieoff.ab=0.1/time.scale) ## incubation time increase, ~pertussis
 State1<-c(nat=1,abt=0)
 Time=seq(0,10*time.scale,by=round(7/(2*time.scale),1))
Pars2<-c(growth.na=0.1,killed.na=1,feed.ab=0.1,dieoff.ab=0.1)/time.scale ##Fast-acting disease, approx btv
State2<-c(nat=20,abt=0)
Pars3<-c(growth.na=1,killed.na=0.5,feed.ab=1.5,dieoff.ab=1)/time.scale ##Chronic infection with acute phase, ~TB
State3<-c(nat=1,abt=0)
Pars4<-c(growth.na=1,killed.na=1.5,feed.ab=2,dieoff.ab=1)  ##Trying to get as similar curves as possible, but skipping this for now....
State4<-c(nat=1,abt=0)



ab.na.diffmod<-function(Time,y=State,parms=Pars){
  with(as.list(c(y,parms)),{
      dna=nat*(growth.na-killed.na*abt)
        dab=feed.ab*nat-dieoff.ab*abt
        return(list(c(dna,dab)))
      })}
diseaseType1=list(parms=Pars1,y=State1)
diseaseType2=list(parms=Pars2,y=State2)
diseaseType3=list(parms=Pars3,y=State3)



 tmp<-ode(func = ab.na.diffmod, y = State1, parms = Pars1, times = Time)
par(mfrow=c(2,1))
 plot(tmp[,1],tmp[,2],type="l")
 plot(tmp[,1],tmp[,3],type="l")


LotkaVolterra.Fun<-function(disease=diseaseType1){
  require(deSolve)
  test.values.lookup<-do.call("ode",c(list(func = ab.na.diffmod,times=seq(0,1000,by=0.5)), disease))
  
  
  function(time){
    time<-ifelse(is.finite(time)&time>max(test.values.lookup[,1]),max(test.values.lookup[,1]),time)
    t1<-approx(test.values.lookup[,1],test.values.lookup[,2],time)$y
    t2<-approx(test.values.lookup[,1],test.values.lookup[,3],time)$y
    t1[!is.finite(time)]<-NA
    t2[!is.finite(time)]<-NA ###This should be modified to exisiting-but-low once we have a latent class model
   return(cbind(t1,t2))
  }
}

par(mfrow=c(3,1))
plot(seq(1,10,by=0.1),LotkaVolterra.Fun(diseaseType1)(seq(1,10,by=0.1))[,1],type="l")
lines(seq(1,10,by=0.1),LotkaVolterra.Fun(diseaseType1)(seq(1,10,by=0.1))[,2],type="l",col="red")

plot(seq(1,20,by=0.1),LotkaVolterra.Fun(diseaseType2)(seq(1,20,by=0.1))[,1],type="l")
lines(seq(1,20,by=0.1),LotkaVolterra.Fun(diseaseType2)(seq(1,20,by=0.1))[,2]*10,type="l",col="red")

plot(seq(1,20,by=0.1),LotkaVolterra.Fun(diseaseType3)(seq(1,20,by=0.1))[,1],type="l",ylim=c(0,3))
lines(seq(1,20,by=0.1),LotkaVolterra.Fun(diseaseType3)(seq(1,20,by=0.1))[,2],type="l",col="red")


