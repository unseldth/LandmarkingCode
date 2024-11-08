
source(paste0(runfiles.path,"AnaLandmark.R"),local = T) 
cum.haz<-readRDS(paste0(runfiles.path,"cumHazSimu.Rds"))


tmat<-transMat(x = list(c(1, 3:5)+1, c(2,3:5)+1,c(1,3:5)+1, c(), c(), c()))
parsL<-data.frame(xi=9,omega = 4.3,alpha = -1)

# landmark parameters
LMps<-seq(1,13)
tt <- seq(min(LMps),max(LMps),by=1)
nt<-length(tt)

##########################################################################################################
## Simulation ############################################################################################
## ######################################################################################################
## Scenario 1 ##############################

simMSM<-function(iteration,nit,nparallel){
  simSeed<-(nparallel-1)*nit+iteration
  print(simSeed)
  set.seed(simSeed)
  out<-mssample(Haz=cum.haz,
                trans=tmat, 
                M=n, clock="forward", 
                history=list(state=sample(c(1,2),n,replace=T,prob=c(0.6,0.4)), time=0),
                output="data",
                do.trace = NULL)
  out<-data.table(as.data.frame(do.call("cbind",out)))
  return(cbind("it"=iteration,out))
}


## Scenario 2 & 3 #########################
appendScenarios<-function(out1){
  ## Scenario 2-----------------------------------------------------------------------------------------
  output<-copy(out1)
  output[,"nr":=paste(it,id,sep="_")]
  LT<-pmax(0,rsn(n, xi=parsL[1,1],omega = parsL[1,2],alpha = parsL[1,3])) # sn= Skew-Normal Distribution
  output2<-output[with(output, order(it, id)),]
  rows.per.subj<-rle(as.character(output$nr))
  dLT<-data.table("nr"=rep(rows.per.subj$values,rows.per.subj$lengths), "LT"=rep(LT, rows.per.subj$lengths))
  temp<-cbind("LT"=dLT$LT, output)
  temp.complete<-temp[LT<Tstart] # completely observed obs
  temp.part<-temp[LT>=Tstart & LT<=Tstop] # only observed from LT on
  temp.part[,"Tstart":=LT]
  
  output2<-rbind(temp.complete,temp.part)
  output2<-output2[with(output2, order(it, id,Tstart)),]
  output2[,c("nr","LT","duration"):=list(NULL,NULL,Tstop-Tstart)]
  
  out2<-cbind(scenario=rep(2,nrow(output2)),output2)
  
  ## Scenario 3-----------------------------------------------------------------------------------------
  simDataLT4<-copy(output2)
  LT<-4 # left censoring at week 4
  simDataLT4<-output2[Tstop>LT] # completed transitions before week 4 are administratively censored
  simDataLT4[,"Tstart":=ifelse(Tstart<=LT & Tstop>LT ,LT,Tstart)] # partly observed transitions
  
  out3<-cbind(scenario=rep(3,nrow(simDataLT4)),simDataLT4)
  
  out1<-cbind(scenario=rep(1,nrow(out1)),out1)
  simData<-rbind(out1,out2,out3)
  
  setnames(simData,old="id",new = "nummer")
  
  return(simData)
}




##########################################################################################################
## Landmark analysis #####################################################################################
## ######################################################################################################
## Preparation ##########################################
tmatCompete <- trans.comprisk(3, names = c("LiveBirth","SAB","ETOP"))
f<-'Surv(Tstart,Tstop, status) ~exposure1.1 +exposure1.2+exposure1.3+ exposure2.1 + exposure2.2 + exposure2.3 + strata(trans)'
#####################################################################################################
## Estimation #####################################################################################
## Crude ----------------------------------------------------------------------------------
lmEstsSeparate<-function(LMlong){
  probs<-NULL
  coefs<-NULL
  for(sc in scenariosSeperate){
    LMlong.sc<-subset(LMlong, scenario==sc)
    LMps<-unique(LMlong.sc$LM) # only choose Landmark points where there is data
    res<-trans.probs.LM(LMlong.sc,f, LMp=LMps)
    coefs<-rbind(coefs,cbind(scenario=sc, res[[1]]))
    probs<-rbind(probs,cbind(scenario=sc, res[[2]]))
  }
  return(list(coefs,probs))
}


## Supermodel ----------------------------------------------------------------------------------
## Formulas##
g1 <- function(t) t-min(LMps)
g2 <- function(t) (t-min(LMps))^2
exposure.terms<-c("exposure1.t1","exposure1.t2","exposure2.t1", "exposure2.t2")
fSuper<-"Surv(Tstart,Tstop, status) ~exposure1.1 +exposure1.2+exposure1.3+ exposure2.1 + exposure2.2+ exposure2.3 + LM1.1+LM1.2+LM1.3 +LM2.1+LM2.2+LM2.3 + strata(trans)"

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
  for(sc in scenarios){
    LMlong.sc<-subset(LMlong, scenario==sc)
    tt<-LMps
    #tt<-unique(LMlong.sc$LM) # only choose Landmark points where there is data
    # just for illustration: supermodel wit all interactions
    LMsupercox.all  <- coxph(as.formula(fSuper.interactions.all), data = LMlong.sc, method = "breslow")
    coxRes.all<-summary(LMsupercox.all)$coefficients
    coxRes.all<-data.table(scenario=sc, variable=rownames(coxRes.all),coxRes.all[,1:4])
    coefsSupermodelAll<-rbind(coefsSupermodelAll,coxRes.all)
    
    ## Proper Model ##############################
    LMsupercox  <- coxph(as.formula(fSuper), data = LMlong.sc, method = "breslow")
    coxRes<-summary(LMsupercox)$coefficients
    coxRes<-data.table(scenario=sc, variable=rownames(coxRes),coxRes[,1:4])
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
    probs<-cbind(scenario=rep(sc,nrow(probs)),probs)
    probsSupermodel<-rbind(probsSupermodel,probs)
  }
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

msEstimates<-function(sc,msData){
  c0<-coxph(Surv(Tstart, Tstop, status)~strata(trans),data=msData[scenario==sc], method="breslow")
  msf0 <- msfit(object = c0, vartype = "greenwood", trans =tmat)
  probslist<-lapply(LMps, mspred, msf0=msf0)
  probs<-rbindlist(probslist)
  return(probs)
}





##############################################################################################
#### create tables for description of landmark datasets #######################################
#############################################################################################
lmCounts<-function(sc,LMlong){
  LMlong2<-data.table(copy(LMlong))
  data<-LMlong2[status==1 & scenario==sc]
  nrisk<-data[,list("Never before"=sum(exposure=="0"),"Currently"=sum(exposure=="1"),"Previously"=sum(exposure=="2" )), by=LM]
  nrisk.m <- melt(nrisk,id.vars = "LM") 
  
  nevent<-data[,list("Life Birth"=sum(trans==1),"SAB"=sum(trans==2), "ETOP"=sum(trans==3)), by=LM]
  nevent.m <- melt(nevent,id.vars = "LM") 
  
  nriskevent<-rbind(cbind(nrisk.m,type="nrisk"),cbind(nevent.m,type="nevent"))
}













