library(data.table)
library(mstate)
library(doParallel)



# Landmark parameters
LMps<-seq(5,13)
horizon<-45


# expand a covariate to transition specific covariates
# inputs: dataset in multi-state form (obs), covariates (covs)
# ouput: expanded dataset
expand.exposure<-function(obs,cov){
  attr(obs, "trans") <- tmatCompete
  class(obs) <- c("msdata", "data.frame")
  cbind(expand.covs(obs, cov, longnames = FALSE), strata=1:nevents)
}


## function for making coefficient tables (output) with CIs from Cox summary (input)
confint.beta<-function(x,cov_values){
  res<-matrix(0,nrow=length(cov_values),ncol=3)
  for(i in 1:length(cov_values)){
    HR <-signif(x$coef[cov_values[i],"exp(coef)"], digits=3);#exp(beta)
    HR.confint.lower <- round(x$conf.int[cov_values[i],"lower .95"], 3)
    HR.confint.upper <- round(x$conf.int[cov_values[i],"upper .95"],3)
    res[i,]<-c(HR, HR.confint.lower,HR.confint.upper)
  }
  colnames(res)<-c("HR","lower .95","upper .95")
  rownames(res)<-cov_values
  return(res)
}


###############################################################
## Functions for data manipulation and prepartion #############
###############################################################

# inouts: data in multistate format (data.MSM), covariates of interes (covs),
# unique number that identifies the subjects in the original dataset (idOrig)
# ouput: stacked landmark data set
prepare.LMData<-function(data.MSM,covs=NULL, idOrig="idChild"){
  covs<-c("exposure",covs)
  data.MSM<-data.table(data.MSM)

  # data.MSM[,"exposure":=ifelse(from==0,"0","1")]
  # 
  # data.MSM[,c("entry","exit","event"):=list(min(Tstart),max(Tstop),max(as.numeric(to))),by=idOrig]
  ############################################################################
  LMall<-foreach(j = 1:length(LMps),.combine=rbind)%do%{
    LM<-LMps[j]

    
    # Select the subjects at risk at the landmark
    LMdata<-data.MSM[entry<=LM & exit>LM]
    
    # Define exposure state at landmark (1 if woman was exposed before at least once at LM)
    LMdata[!is.na(from_1),"exposure":=ifelse(from_1<=LM,"1","0")]
    LMdata[is.na(from_1),"exposure":="0"]
    table(LMdata$exposure)
    LMdata[exposure=="1"]
    LMdata[exposure=="0"]
    
    
    
    # Define event status
    LMdata[,"Tstop":=exit]
    LMdata[,"status":=fcase(to==3,"1",
                           to==4,"2",
                           to==5,"3",
                           to==6,"4",
                           default=NA)]
    # Consider events that occur within the prediction horizon (other women are considered event-free at horizon)
    LMdata[exit>horizon,c("Tstop","status"):=list(horizon,NA)]
    
    LMdata[,c("stat1","stat2","stat3","stat4"):=list(as.numeric(status=="1"),
                                             as.numeric(status=="2"),
                                             as.numeric(status=="3"),
                                             as.numeric(status=="4"))]
    head(LMdata)

    LMdata[,"LM":=LM]
    LMlong <- msprep(time =c(NA,rep("Tstop",nevents)), status =c(NA,paste0("stat",1:nevents)), data = LMdata,  
                     keep = c("idWoman", "LM","entry",covs),
                     trans = tmatCompete,
                     id=idOrig)
    LMlong
    return(LMlong)
  }
  

  rownames(LMall)<-NULL
  attr(LMall, "trans") <- tmatCompete

  LMall$Tstart<-LMall$LM # The entry/left truncation time is the landmark time
  LMall$time<-NULL

  LMlongcovs <- expand.covs(LMall, covs,longnames = F)
  head(LMlongcovs)
  return(LMlongcovs)
}


########################################################################################
## complementary log-log transformation of confidence interval Pgj(LM,t) ###############
########################################################################################
# input: matrix of probability estimates and standard errors (pts),
# vectors of covariates, landmarks, and outcome states of interest (exposure, 
# LM, pstats)
# output: predictions with log-log confidence intervals
pt.pred<-function(exposure,pts,LM,pstats=(2+seq(0,nevents-1))){
  pt<-pts[pts$exposure==exposure,]
  ps<-matrix(as.numeric(cbind(pt[,paste0("pstate",pstats)], pt[,paste0("se",pstats)])),ncol=2)
  pt.ci<-data.frame("LM"=rep(LM,length(pstats)),
                   "exposure"=rep(exposure,length(pstats)),
                   "trans"=1:length(pstats),
                   signif(do.call("rbind",apply(ps,1,confint.trans)),4)
             )
  return(pt.ci)
}

confint.trans <- function(p.se) {
  p<-p.se[1]
  se<-p.se[2]
  alpha <- qnorm(0.975)
  lower <- 1 - (1 - p)^(exp(alpha * (se/((1 - p) * log(1 - p)))))
  upper <- 1 - (1 - p)^(exp(-alpha * (se/((1 - p) * log(1 - p)))))
  return(data.frame(CIF=p,lower, upper))
}


##############################################################
# calculate transition probabilities from msfit objects ######
##############################################################

# input: (Cox regression) formula (transf), vector of covariates in expandend form (cov_values)
trans.probs.LM2<-function(transf,cov_values=c("exposure.1","exposure.2","exposure.3")){
  probs<-NULL
  ci.betas<-NULL
  for(landmark in LMps){
    s1<- subset(LMlong,LM==landmark)
    c1 <- coxph(as.formula(transf), data = subset(LMlong,LM==landmark), method = "breslow")

    ci.beta<-confint.beta(x=summary(c1),cov_values=cov_values) 
    ci.betas<-rbind(ci.betas,
                    data.table("LM"=rep(landmark,nrow(ci.beta)),"Variable"=rownames(ci.beta),ci.beta))
    ci.betas
    
    msf.never <- msfit(c1, nexposure, trans = tmatCompete)
    msf.expo <- msfit(c1, before.expo, trans = tmatCompete)
    
    pts<-cbind(exposure=0:1,
               rbind(tail(probtrans(msf.never, predt=landmark)[[1]],n=1),
                     tail(probtrans(msf.expo, predt=landmark)[[1]],n=1)))
    
    probs<-rbind(probs,data.table(do.call("rbind",lapply(0:1, pt.pred,pts=pts,LM=landmark))))
  }
  return(list(ci.betas,probs))
}


