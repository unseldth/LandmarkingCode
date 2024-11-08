## Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"HelpFunctionsAnaGyrase.R"), local = T)


theme_set(theme_bw(base_size=16))

LMlong<-readRDS(paste0(data.path,"/LMlong.Rds"))


#-------------------------------------------------------------------------------------------

############################################################################################
## Landmark supermodel #####################################################################
############################################################################################

# Time points at which we want predictions
tt <- seq(min(LMps),max(LMps),by=1)
nt<-length(tt)

## Model #####################################################################################
g1 <- function(t) t-min(LMps)
g2 <- function(t) (t-min(LMps))^2
LMlong$LM1 <- g1(LMlong$LM)
LMlong$LM2 <- g2(LMlong$LM)


LMlong$exposure1.t1<-as.numeric(LMlong$exposure=="1")* g1(LMlong$LM) #interaction with s
LMlong$exposure1.t2<-as.numeric(LMlong$exposure=="1")* g2(LMlong$LM) #interaction with s^2

attr(LMlong, "trans") <- tmatCompete
LMlong <- expand.covs(LMlong, c("LM1","LM2", "exposure1.t1","exposure1.t2"))
LMlong[,c("LM1","LM2", "exposure1.t1","exposure1.t2")]<-NULL


# paste(paste0("exposure.",1:nevents),"+",paste0("LM1.",1:nevents),paste0("LM2.",1:nevents),collapse="+")
f<-paste0('Surv(Tstart,Tstop, status) ~exposure.1 +exposure.2+exposure.3+LM1.1 + LM1.2 + LM1.3+',
          'LM2.1+LM2.2+LM2.3+ strata(trans)+cluster(idChild)')
LMsupercox  <- coxph(as.formula(f),data = LMlong, method = "breslow", robust=T)
summary(LMsupercox)

# without LM^2
f<-paste0('Surv(Tstart,Tstop, status) ~exposure.1 +exposure.2+exposure.3+LM1.1 + LM1.2 + LM1.3+',
          'strata(trans)+cluster(idChild)')
LMsupercox1  <- coxph(as.formula(f),data = LMlong, method = "breslow", robust=T)
summary(LMsupercox1)

tabCox<-summary(LMsupercox1)
resCox<-confint.beta(tabCox,cov_values=rownames(tabCox$coefficients))
temp<-apply(resCox,1, function(i) paste0(i[1]," (",i[2],",",i[3],")"))
temp
#"Age"=rep(NA,3),"moderate"=rep(NA,3)
tabprintadj<-t(data.table("Fluoroquinolone"=temp[1:3],
                          "$\\gamma_1$"=temp[4:6]))
colnames(tabprintadj)<-event.lbls[-4]


## with interaction terms #################################
# quadratic and linear
exposure.interactions<-paste0("exposure1.t1.1+exposure1.t1.2+exposure1.t1.3+exposure1.t2.1+exposure1.t2.2+",
                              "exposure1.t2.3")
f2<-paste(f,"+",exposure.interactions)
LMsupercox2  <- coxph(as.formula(f2), data = LMlong, method = "breslow", robust=T)
summary(LMsupercox2)

# only linear
exposure.interactions.linear<-"exposure1.t1.1+exposure1.t1.2+exposure1.t1.3"
f2<-paste(f,"+",exposure.interactions.linear)
LMsupercox2  <- coxph(as.formula(f2), data = LMlong, method = "breslow", robust=T)
summary(LMsupercox2)
cbind(summary(LMsupercox2)$coefficients[,c("exp(coef)","Pr(>|z|)")],
      summary(LMsupercox2)$conf.int[,c("lower .95","upper .95")])
# exposure interactions are not significant


###########################################################################
## Unadjusted model (with just fluoroquinolone exposure) ###################
###########################################################################


# LM^2 not significant but is kept because it gives smoother curves and it is needed for the adjusted model
# estimated coefficients #####################################################
f<-paste0('Surv(Tstart,Tstop, status) ~exposure.1 +exposure.2+exposure.3+LM1.1 + LM1.2 + LM1.3+',
          'LM2.1+LM2.2+LM2.3+ strata(trans)+cluster(idChild)')
LMsupercox  <- coxph(as.formula(f),data = LMlong, method = "breslow", robust=T)

printCoefsCox<-function(summaryCox){
  resCox<- confint.beta(summaryCox,cov_values=rownames(summaryCox$coefficients))
  temp<-apply(resCox,1, function(i) paste0(i[1]," (",i[2],",",i[3],")"))
  tabprint<-t(data.table("Fluoroquinolone"=temp[1:3], "$\\gamma_1$"=temp[4:6],"$\\gamma_2$"=temp[7:9]))
  colnames(tabprint)<-event.lbls[-4]
  temp<-xtable(tabprint, caption=paste("Estimated regression parameters of the landmark supermodel",
                                       "with both quadratic and linear landmark main effects on an exponential scale. $95\\%$ Confidence Interval in brackets"), 
               label="tab:supermodel coefficients")
  print(temp,sanitize.text.function=function(x){x})
}
printCoefsCox(summary(LMsupercox))



## Prediction ################################################################################
Fw.never<-matrix(NA,nrow = nt, ncol=1+2*(nevents+1))
Fw.expo<-matrix(NA,nrow = nt, ncol=1+2*(nevents+1))
probsSuper<-NULL

for(i in 1:nt){
  never<-data.frame(exposure = factor(rep(0, nevents),levels=0:1), LM1=rep(g1(tt[i]),nevents),LM2=rep(g2(tt[i]),nevents), 
                    trans=1:nevents)
  expo<-data.frame(exposure = factor(rep(1, nevents),levels=0:1), LM1=rep(g1(tt[i]),nevents),LM2=rep(g2(tt[i]),nevents),
                   trans=1:nevents)
  
  never<-expand.exposure(never, cov=c("exposure","LM1","LM2"))
  expo<-expand.exposure(expo, cov=c("exposure","LM1","LM2"))
  
  msf.never<- msfit(LMsupercox, never, trans = tmatCompete)
  msf.expo<- msfit(LMsupercox, expo, trans = tmatCompete)
  
  
  Fw.never[i,]<-as.numeric(tail(probtrans(msf.never, predt=tt[i])[[1]],n=1))
  Fw.expo[i,]<-as.numeric(tail(probtrans(msf.expo, predt=tt[i])[[1]],n=1))
  
  pts<-cbind(exposure=0:1,
             rbind(tail(probtrans(msf.never, predt=tt[i])[[1]],n=1),
                   tail(probtrans(msf.expo, predt=tt[i])[[1]],n=1)))
  
  probsSuper<-rbind(probsSuper,data.table(do.call("rbind",lapply(0:1, pt.pred,pts=pts,LM=tt[i]))))
  
}
probsSuper$exposure<-factor(probsSuper$exposure, levels=0:1)

saveRDS(probsSuper,file = paste0(wdresults,"probsSuper.Rds"))







######################################################
########## Adujsted for age ##########################
######################################################

LMlong$age.gr1.t1<-as.numeric(LMlong$age.gr=="<=20 years")* g1(LMlong$LM) #interaction with s
LMlong$age.gr1.t2<-as.numeric(LMlong$age.gr=="<=20 years")* g2(LMlong$LM) #interaction with s^2
LMlong$age.gr2.t1<-as.numeric(LMlong$age.gr==">35 years")* g1(LMlong$LM) #interaction with s
LMlong$age.gr2.t2<-as.numeric(LMlong$age.gr==">35 years")* g2(LMlong$LM) #interaction with s^2

LMlong2 <- expand.covs(LMlong, c("age.gr1.t1","age.gr1.t2", "age.gr2.t1","age.gr2.t2"))

age.constant<-"age.gr1.1+age.gr1.2+age.gr1.3+age.gr2.1+age.gr2.2+age.gr2.3" # constant terms to be kept in the model
age.interactions<-paste0("age.gr1.t1.1+age.gr1.t1.2+age.gr1.t1.3+age.gr1.t2.1+age.gr1.t2.2+age.gr1.t2.3+age.gr2.t1.1+",
                         "age.gr2.t1.2+age.gr2.t1.3+age.gr2.t2.1+age.gr2.t2.2+age.gr2.t2.3")



#### backward selection of interaction terms based on Wald tests ############################################
f<-paste0('Surv(Tstart,Tstop, status) ~exposure.1 +exposure.2+exposure.3+LM1.1 + LM1.2 + LM1.3 +LM2.1+LM2.2+LM2.3 + strata(trans)+',
          'cluster(idChild)')
f2<-paste0(f,"+",age.constant) # formula part that should be kept anyway
interactions<-paste0(exposure.interactions,"+",age.interactions)
interactions.names<-unlist(strsplit(interactions,split="\\+"))

pval<-1
while(pval>0.05){
  fadj<-paste0(f2,"+",paste(interactions.names,collapse="+"))
  LMsupercoxadj1<- coxph(as.formula(fadj),data = LMlong2, method = "breslow", robust=T)
  coefs<-summary(LMsupercoxadj1)$coefficients
  pval<-max(coefs[interactions.names,"Pr(>|z|)"])
  pval
  if(pval<=0.05) return()
  to.remove<-names(which.max(coefs[interactions.names,"Pr(>|z|)"]))
  print(c(pval,to.remove))
  interactions.names<-setdiff(interactions.names,to.remove)
  cat("Interactions:", interactions.names)
}


interactions.names

# there remains 1 interactions: 
# age.gr2.t2.3 (interaction effect of the exposure with quadratic term of landmark time on ETOP),
# so we continue with this (plus the corresponding linear term)

fadj<-paste0("Surv(Tstart,Tstop, status) ~exposure.1 +exposure.2+exposure.3+LM1.1 + LM1.2 + LM1.3 +LM2.1+LM2.2+LM2.3 +",
             "age.gr1.1 + age.gr1.2 + age.gr1.3 + age.gr2.1 + age.gr2.2 + age.gr2.3 +",
             "age.gr2.t1.3 + age.gr2.t2.3 + strata(trans) + cluster(idChild)")
LMsupercoxadj1<- coxph(as.formula(fadj),data = LMlong2, method = "breslow", robust=T)
coefs<-summary(LMsupercoxadj1)$coefficients

####################################################################
### Table 3 Hazard ratios in the landmark supermodel ###############
####################################################################

tabCox<-summary(LMsupercoxadj1)
tabCox
resCox<-data.table(confint.beta(tabCox,cov_values=setdiff(rownames(tabCox$coefficients),"exposure")))
resCox

temp<-with(resCox,sprintf("%.3f (%.3f,%.3f)", HR,`lower .95`,`upper .95`))
tabprintadj<-t(data.table("Fluoroquinolone"=temp[1:3],
                          "Age $\\leq 20$ years "=temp[10:12], "Age $>35$ years"=temp[13:15],
                          "Age $>35$ years $\\cdot$ s"=c(rep(NA,2),temp[16]),
                          "$\\exp(\\gamma_1^{(j)})$"=temp[4:6],"$\\exp(\\gamma_2^{(j)})$"=temp[7:9]))
colnames(tabprintadj)<-event.lbls[-4]
xtab_supermodel_interactions<-xtable(tabprintadj,digits=4, caption="Estimated regression parameters of the landmark supermodel. $95\\%$ Confidence Interval in brackets", 
                                     label="tab:supermodel coefficientsadj")
print(xtab_supermodel_interactions,sanitize.text.function=function(x){x})



## for a comparison: with only linear landmark effects --------------------------------
fadj<-paste0("Surv(Tstart,Tstop, status) ~exposure.1 +exposure.2+exposure.3+LM1.1 + LM1.2 + LM1.3+",
             "age.gr1.1+age.gr1.2+age.gr1.3+age.gr2.1+age.gr2.2+age.gr2.3+age.gr2.t1.3+strata(trans)+cluster(idChild)")
LMsupercoxadj2<- coxph(as.formula(fadj),data = LMlong2, method = "breslow", robust=T)
summary(LMsupercoxadj2)
coefs

## models of significant interaction of covariate age for higher age group ---------
coefs["age.gr2.3","coef"]
coefs["age.gr2.t1.3","coef"]
gam.age2.3<-function(s) exp(coefs["age.gr2.3","coef"]+coefs["age.gr2.t1.3","coef"]*g1(s))
gam.age2.3(min(LMps))
gam.age2.3(6)
gam.age2.3(max(LMps))
par(mar=c(5,4,4,2) + 1)
curve(gam.age2.3(x),min(LMps),max(LMps),xlab="Landmark",ylab=TeX("exp($\\hat{\\beta}_{ageHigh}^{(5)}(s)$)"))
abline(h=0)

## Prediction ######################
#agr=age.lbls[1];i=1
lmSuperAge<-function(agr){
  probsSuperadj<-NULL
  Fw.never<-matrix(NA,nrow = nt, ncol=1+2*(nevents+1))
  Fw.expo<-matrix(NA,nrow =  nt, ncol=1+2*(nevents+1))
  for(i in 1:nt){
    never<-data.frame(exposure = factor(rep(0, nevents),levels=0:1), 
                      LM1=rep(g1(tt[i]),nevents), LM2=rep(g2(tt[i]),nevents), 
                      age.gr=factor(rep(agr, nevents),levels=age.lbls),
                      age.gr2.t1=as.numeric(agr==">35 years")*g1(tt[i]),
                      trans=1:nevents)
    expo<-data.frame(exposure = factor(rep(1, nevents),levels=0:1),
                     LM1=rep(g1(tt[i]),nevents),LM2=rep(g2(tt[i]),nevents), 
                     age.gr=factor(rep(agr, nevents),levels=age.lbls),
                     age.gr2.t1=as.numeric(agr==">35 years")*g1(tt[i]),
                     trans=1:nevents)
    
    
    never<-expand.exposure(never, cov=c("exposure","LM1","LM2","age.gr","age.gr2.t1"))
    expo<-expand.exposure(expo, cov=c("exposure","LM1","LM2","age.gr","age.gr2.t1"))
    
    msf.never<- msfit(LMsupercoxadj1, never, trans = tmatCompete)
    msf.expo<- msfit(LMsupercoxadj1, expo, trans = tmatCompete)
    
    Fw.never[i,]<-as.numeric(tail(probtrans(msf.never, predt=tt[i])[[1]],n=1))
    Fw.expo[i,]<-as.numeric(tail(probtrans(msf.expo, predt=tt[i])[[1]],n=1))
    
    pts<-cbind(exposure=0:1,
               rbind(tail(probtrans(msf.never, predt=tt[i])[[1]],n=1),
                     tail(probtrans(msf.expo, predt=tt[i])[[1]],n=1)))
    
    probsSuperadj<-rbind(probsSuperadj,data.table(do.call("rbind",lapply(0:1, pt.pred,pts=pts,
                                                                         pstats=(2+seq(0,nevents-1)),
                                                                         LM=tt[i]))))
    
  }
  
  probsSuperadj$exposure<-factor(probsSuperadj$exposure, levels=0:1)
  return(probsSuperadj)
}

fadj<-paste0(f2, "+age.gr2.t1.3")
LMsupercoxadj1<- coxph(as.formula(fadj),data = LMlong2, method = "breslow", robust=T)

probsAgelist<-lapply(age.lbls,lmSuperAge)
probsAge<-rbindlist(probsAgelist,idcol="Age")
probsAge[,"exposure":=factor(exposure,levels=c("0","1"),labels=c("No prior exposure","Prior exposure"))]
probsAge[,"Age":=factor(Age,levels=c(2,1,3),labels=age.lbls[c(2,1,3)])]

saveRDS(probsAge,file=paste0(wdresults,"probsSuperAdj.Rds"))



##### Smooth function of landmark main effects ###################
coefsLMsupercox<-summary(LMsupercoxadj1)$coefficients
coefsLMsupercox
gam1<-function(s) coefsLMsupercox[paste0("LM1.",1),"coef"]*g1(s)+ coefsLMsupercox[paste0("LM2.",1),"coef"]*g2(s)
gam2<-function(s) coefsLMsupercox[paste0("LM1.",2),"coef"]*g1(s)+ coefsLMsupercox[paste0("LM2.",2),"coef"]*g2(s)
gam3<-function(s) coefsLMsupercox[paste0("LM1.",3),"coef"]*g1(s)+ coefsLMsupercox[paste0("LM2.",3),"coef"]*g2(s)
curve(exp(gam1(x)),min(LMps),max(LMps),ylim=c(0.75,1.2),xlab="Landmark",ylab=TeX("exp($\\hat{\\gamma}_j(s)$)"))
curve(exp(gam2(x)),min(LMps),max(LMps),add=T,lty=2, col="dodgerblue4")
curve(exp(gam3(x)),min(LMps),max(LMps),add=T,lty=3, lwd=1.5, col="darkolivegreen4")
legend("topleft",legend=c("Live Birth","SAB","ETOP"),lty=1:3,bty="n",cex = 0.8, 
       col = c("black","dodgerblue4","darkolivegreen4"))




