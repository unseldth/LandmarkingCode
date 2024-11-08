source(paste0(runfiles.path,"AnaLandmarkScs.R"),local = T) 

## simulation parameters
cum.hazScs<-readRDS(paste0(runfiles.path,"cumHazSimuScs.Rds"))
tmat<-transMat(x = list(c(1, 3:5)+1, c(2,3:5)+1,c(1,3:5)+1, c(), c(), c()))
tsim<-seq(0.00,50,0.05)
cutt<-45 # no events after this week
cuttSAB<-23
LTpars=list(list(distr="SN",xi=9,omega = 4.3,alpha = -1),
            list(distr="Exp",rate=0.1),
            list(distr="Exp",rate=0.5),
            list(distr="Exp",rate=1.5))
            

weibull.cumhaz<-function(t,shape,scale){scale*t^shape}
pars.SAB<-data.frame(shape=c(0.45),scale=c(0.03))
HexpETOP<-function(t,steps,t.limits) {
  steps[1]*(min(t,t.limits[2])-t.limits[1])*1*(t>=t.limits[1])+
    steps[2]*(min(t,t.limits[3])-t.limits[2])*1*(t>t.limits[2])+
    steps[3]*(min(t,t.limits[4])-t.limits[3])*1*(t>t.limits[3])+
    steps[4]*(min(t,t.limits[5])-t.limits[4])*1*(t>t.limits[4])+
    steps[5]*(min(t,t.limits[6])-t.limits[5])*1*(t>t.limits[5])+
    steps[6]*(min(t,t.limits[7])-t.limits[6])*1*(t>t.limits[6])
}
tm <- c(0,2,5,12,23,30,cutt)
ta <- c(1/1000,1/100,4/100,1/1000,1/10000,0)/2
pars.ETOP<-matrix(ta,nrow=1,byrow=T)
myHexp<-function(x,rate){ rate*x}


# landmark parameters
LMps<-seq(1,13)
tt <- seq(min(LMps),max(LMps),by=1)
nt<-length(tt)

##########################################################################################################
## Parameters ############################################################################################
## ######################################################################################################

specifyCumHaz<-function(ParsSetting){
  ##SAB
  A04<-sapply(pmin(tsim,cuttSAB),weibull.cumhaz, shape=pars.SAB[1,1],scale=pars.SAB[1,2])
  A14<-ParsSetting$HRSABCur*A04
  A24<-ParsSetting$HRSABPrev*A04
  
  ##ETOP
  A05<-sapply(pmin(tsim,cutt),  HexpETOP,steps=pars.ETOP[1,],t.limits=tm)
  A15<-ParsSetting$HRETOPCur*A05
  A25<-ParsSetting$HRETOPPrev*A05
  
  ##Cur->Prev
  A12<-sapply(pmin(tsim,cutt),myHexp, rate=ParsSetting$RateCurPrev)
  
  cum.haz<-rbind(cum.hazScs,
                 data.frame(time=tsim,Haz=A04,trans=as.integer(3)),
                 data.frame(time=tsim,Haz=A05,trans=as.integer(4)),
                 data.frame(time=tsim,Haz=A12,trans=as.integer(5)),
                 data.frame(time=tsim,Haz=A14,trans=as.integer(7)),
                 data.frame(time=tsim,Haz=A15,trans=as.integer(8)),
                 data.frame(time=tsim,Haz=A24,trans=as.integer(11)),
                 data.frame(time=tsim,Haz=A25,trans=as.integer(12))
  )
  
  
}


##########################################################################################################
## Simulation ############################################################################################
## ######################################################################################################
## Scenario 1 ##############################

simMSM<-function(iteration,ParsSetting,nit,nparallel){
  cum.haz<-specifyCumHaz(ParsSetting)
  simSeed<-(nparallel-1)*nit+iteration
  print(simSeed)
  set.seed(simSeed)
  out<-mssample(Haz=cum.haz,
                trans=tmat, 
                M=n, clock="forward", 
                history=list(state=sample(c(1,2),n,replace=T,prob=c(ParsSetting$pStartNever,1-ParsSetting$pStartNever)), 
                             time=0),
                output="data",
                do.trace = NULL)
  out<-data.table(as.data.frame(do.call("cbind",out)))
  out1<-cbind("it"=iteration,out)
  out1[Tstop==Inf,c("Tstop","status"):=list(45,0)]
  out1[,"duration":=Tstop-Tstart]

  return(out1)
}


## Scenario 2 & 3 #########################
appendScenarios<-function(out1,LTi, iteration,nit,nparallel){
  ## Scenario 2-----------------------------------------------------------------------------------------
  output<-copy(out1)

  set.seed((nparallel-1)*nit+iteration)
  LeftTruncDist<-LTpars[[LTi]]
  if(LeftTruncDist$distr=="Exp"){
    LT<-rexp(n,rate=LeftTruncDist$rate)
  }else if(LeftTruncDist$distr=="SN"){
    LT<-pmax(0,rsn(n, xi=LeftTruncDist$xi,omega = LeftTruncDist$omega,alpha = LeftTruncDist$alpha))
  }
  output2<-output[with(output, order(id)),]
  rows.per.subj<-rle(as.character(output$id))
  dLT<-data.table("LT"=rep(LT, rows.per.subj$lengths))
  temp<-cbind("LT"=dLT$LT, output)
  temp.complete<-temp[LT<Tstart] # completely observed obs
  temp.part<-temp[LT>=Tstart & LT<=Tstop] # only observed from LT on
  temp.part[,"Tstart":=LT]
  
  output2<-rbind(temp.complete,temp.part)
  output2<-output2[with(output2, order(id,Tstart)),]
  output2[,"nummer":=id]
  output2[,c("id","LT","duration"):=list(NULL,NULL,Tstop-Tstart)]
  
  return(output2)
}


##########################################################################################################
## Time-dependent covariate ##############################################################################
## ######################################################################################################
## Preparation ##########################################
estsTdc<-function(outcome,TdcData){
  tab_tdc=summary(coxph(Surv(Tstart,Tstop,to==outcome)~factor(from),data=TdcData))$coefficients[,"exp(coef)"]
  return(data.table(exposure=1:2,HR=tab_tdc))
}

##########################################################################################################
## Landmark analysis #####################################################################################
## ######################################################################################################
## Preparation ##########################################
tmatCompete <- trans.comprisk(3, names = c("LiveBirth","SAB","ETOP"))
f<-'Surv(Tstart,Tstop, status) ~exposure1.1 +exposure1.2+exposure1.3+ exposure2.1 + exposure2.2 + exposure2.3 + strata(trans)'
#####################################################################################################
## Estimation #####################################################################################
## Separate ----------------------------------------------------------------------------------
lmEstsSeparate<-function(LMlong){
  LMps<-unique(LMlong$LM) # only choose Landmark points where there is data
  res<-trans.probs.LM(LMlong,f, LMp=LMps)
  coefs<-res[[1]]
  probs<-res[[2]]
  return(list(coefs,probs))
}


## Supermodel ----------------------------------------------------------------------------------
## Formulas##
g1 <- function(t) t-min(LMps)
g2 <- function(t) (t-min(LMps))^2
exposure.terms<-c("exposure1.t1","exposure1.t2","exposure2.t1", "exposure2.t2")
fSuper<-"Surv(Tstart,Tstop, status) ~exposure1.1 +exposure1.2+exposure1.3+ exposure2.1 + exposure2.2+ exposure2.3 + LM1.1+LM1.2+LM1.3 +LM2.1+LM2.2+LM2.3 + strata(trans)"

# with interactions of exposures and landmark times
exposure.interactions.all<-"exposure1.t1.1 +exposure1.t1.2 +exposure1.t1.3 +exposure1.t2.1 +exposure1.t2.2 +exposure1.t2.3 +exposure2.t1.1 +exposure2.t1.2 +exposure2.t1.3 +exposure2.t2.1 +exposure2.t2.2 +exposure2.t2.3"
fSuper.interactions.all<-paste(fSuper,"+",exposure.interactions.all)


lmEstsSupermodel<-function(LMlong){
  
  LMlong$LM1 <- g1(LMlong$LM)
  LMlong$LM2 <- g2(LMlong$LM)
  
  LMlong$exposure1.t1<-as.numeric(LMlong$exposure==1)* g1(LMlong$LM) #interaction with s
  LMlong$exposure1.t2<-as.numeric(LMlong$exposure==1)* g2(LMlong$LM) #interaction with s^2
  
  LMlong$exposure2.t1<-as.numeric(LMlong$exposure==2)* g1(LMlong$LM) #interaction with s
  LMlong$exposure2.t2<-as.numeric(LMlong$exposure==2)* g2(LMlong$LM) #interaction with s^2
  
  LMlong <- expand.covs(LMlong, c("LM1","LM2",exposure.terms))
  
  probsSupermodel<-NULL
  coefsSupermodel<-NULL
  coefsSupermodelAll<-NULL
    tt<-unique(LMlong$LM) # only choose Landmark points where there is data
    # just for illustration: supermodel wit all interactions
    LMsupercox.all  <- coxph(as.formula(fSuper.interactions.all), data = LMlong, method = "breslow")
    coxRes.all<-summary(LMsupercox.all)$coefficients
    coxRes.all<-data.table(variable=rownames(coxRes.all),HR=coxRes.all[,"exp(coef)"])
    coefsSupermodelAll<-rbind(coefsSupermodelAll,coxRes.all)
    
    ## Proper Model ##############################
    LMsupercox  <- coxph(as.formula(fSuper), data = LMlong, method = "breslow")
    coxRes<-summary(LMsupercox)$coefficients
    coxRes<-data.table(variable=rownames(coxRes),HR=coxRes[,"exp(coef)"])
    coefsSupermodel<-rbind(coefsSupermodel,coxRes)
    
    
    ## Prediction ##################################
    Fw.never<-matrix(NA,nrow = nt, ncol=9)
    Fw.cur<-matrix(NA,nrow = nt, ncol=9)
    Fw.prev<-matrix(NA,nrow = nt, ncol=9)
    
    probs<-NULL
    for(i in seq_along(tt)){
      never<-data.frame(exposure = factor(rep(0, 3),levels=0:2), 
                        exposure1.t1=0*g1(tt[i]),
                        exposure1.t2=0*g2(tt[i]),
                        exposure2.t1=0*g1(tt[i]),
                        exposure2.t2=0*g2(tt[i]),
                        LM1=rep(tt[i],3), LM2=rep(g2(tt[i]),3), 
                        trans=1:3)
      cur<-data.frame(exposure = factor(rep(1, 3),levels=0:2),
                      exposure1.t1=1*g1(tt[i]),
                      exposure1.t2=1*g2(tt[i]),
                      exposure2.t1=0*g1(tt[i]),
                      exposure2.t2=0*g2(tt[i]),
                      LM1=rep(tt[i],3), LM2=rep(g2(tt[i]),3), trans=1:3)
      prev<-data.frame(exposure = factor(rep(2, 3),levels=0:2),
                       exposure1.t1=0*g1(tt[i]),
                       exposure1.t2=0*g2(tt[i]),
                       exposure2.t1=1*g1(tt[i]),
                       exposure2.t2=1*g2(tt[i]),
                       LM1=rep(tt[i],3), LM2=rep(g2(tt[i]),3), trans=1:3)
      
      never<-expand.exposure(never, cov=c("LM1","LM2","exposure",exposure.terms))
      cur<-expand.exposure(cur, cov=c("LM1","LM2","exposure",exposure.terms))
      prev<-expand.exposure(prev,cov=c("LM1","LM2","exposure",exposure.terms))
      
      msf.never<- msfit(LMsupercox, never, trans = tmatCompete)
      msf.cur<- msfit(LMsupercox, cur, trans = tmatCompete)
      msf.prev<- msfit(LMsupercox, prev, trans = tmatCompete)
      
      
      pts<-cbind(exposure=0:2,
                 rbind(tail(probtrans(msf.never, predt=tt[i], variance = F)[[1]],n=1),
                       tail(probtrans(msf.cur, predt=tt[i], variance = F)[[1]],n=1),
                       tail(probtrans(msf.prev, predt=tt[i], variance = F)[[1]],n=1)))
      
      probs<-rbind(probs,data.table(do.call("rbind",lapply(0:2, pt.pred,pts=pts,LM=tt[i]))))
      
    }
    probs$exposure<-factor(probs$exposure, levels=0:2)
    probsSupermodel<-rbind(probsSupermodel,probs)

  return(list(coefsSupermodelAll,coefsSupermodel,probsSupermodel))
}


#############################################################################################
### Non-parametric estimates ##################################################################
#############################################################################################
mspred<-function(tti,msf0){
  pt0 <- probtrans(msf0, direction="forward", predt = tti, method = "greenwood",variance = F)
  pts<-cbind(exposure=0:2,
             rbind(tail(pt0[[1]],n=1), # from never exposed before
                   tail(pt0[[2]],n=1), # from currently exposed
                   tail(pt0[[3]],n=1))) # from previouxly exposed
  probs<-data.table(do.call("rbind",lapply(0:2, pt.pred,pts=pts,LM=tti, pstats=4:6)))
  return(probs)
} 

msEstimates<-function(msData){
  c0<-coxph(Surv(Tstart, Tstop, status)~strata(trans),data=msData, method="breslow")
  msf0 <- msfit(object = c0, vartype = "greenwood", trans =tmat)
  probslist<-lapply(LMps, mspred, msf0=msf0)
  probs<-rbindlist(probslist)
  return(probs)
}




##############################################################################################
#### create tables for description of landmark datasets #######################################
#############################################################################################
lmCounts<-function(LMlong){
  LMlong2<-data.table(copy(LMlong))
  data<-LMlong2[status==1]
  nrisk<-data[,list("Never before"=sum(exposure=="0"),"Currently"=sum(exposure=="1"),"Previously"=sum(exposure=="2" )), by=LM]
  nrisk.m <- melt(nrisk,id.vars = "LM") 
  
  nevent<-data[,list("Life Birth"=sum(trans==1),"SAB"=sum(trans==2), "ETOP"=sum(trans==3)), by=LM]
  nevent.m <- melt(nevent,id.vars = "LM") 
  
  nriskevent<-rbind(cbind(nrisk.m,type="nrisk"),cbind(nevent.m,type="nevent"))
}













