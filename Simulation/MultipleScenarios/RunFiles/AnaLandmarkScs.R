library(data.table)
library(mstate)

states<-1:6
state.lbls<-c("NeverExposed","Exposed","NonExposed","Life Birth","SAB","ETOP")
event.lbls<-c("Life birth","SAB","ETOP")
names(states)<-paste("s",state.lbls,sep=".")

tmatCompete <- trans.comprisk(3, names = c("LiveBirth","SAB","ETOP"))

# Landmark parameters
horizon<-45


# expand covariate to transition specific covariates
expand.exposure<-function(obs,cov){
  attr(obs, "trans") <- tmatCompete
  class(obs) <- c("msdata", "data.frame")
  cbind(expand.covs(obs, cov, longnames = FALSE), strata=1:3)
}

## data frames for predicition
never<-data.frame(exposure = factor(rep(0, 3),levels=0:2), trans=1:3)
never<-expand.exposure(never, cov="exposure")

cur<-data.frame(exposure = factor(rep(1, 3),levels=0:2), trans=1:3)
cur<-expand.exposure(cur, cov="exposure")

prev<-data.frame(exposure = factor(rep(2, 3),levels=0:2), trans=1:3)
prev<-expand.exposure(prev, cov="exposure")

####################################################################################
####################################################################################
prepare.LMData<-function(data.MSM,LMp=LMps,base.covs=NULL, idOrig="nummer"){
  covs<-c("exposure",base.covs)
  data.MSM<-data.table(data.MSM)
  ## Data Management ########################
  # Define Exposure Variables
  data.MSM[,"exposure":=ifelse(from==states["s.NeverExposed"],"0",ifelse(from==states["s.Exposed"],"1","2"))]
  
  data.MSM[,c("entry","exit","event"):=list(min(Tstart),max(Tstop),max(as.numeric(to))),by=c(idOrig)]
  ############################################################################
  LMall<-NULL
  for(j in 1:length(LMp)) {
    #j=1
    LM<-LMp[j]
    # Select the subjects at risk at the landmark
    LMdata<-data.MSM[Tstart<=LM & Tstop>LM]
    LMdata[,c("Tstop","event"):=list(ifelse(exit>horizon,horizon,exit),
                                     ifelse(exit>horizon,-1,event)) ]
    
    LMall<-rbind(LMall,cbind(LM=rep(LM, nrow(LMdata)),LMdata))
  }
  
  LMall[,"status":=ifelse(event==6,"3",ifelse(event==5,"2",ifelse(event==4,"1","0")))]
  LMall[,c("stat1","stat2","stat3"):=list(as.numeric(LMall$status=="1"), 
                                          as.numeric(LMall$status=="2"),
                                          as.numeric(LMall$status=="3"))]
  LMall[,"id":=idOrig]
  LMlong <- msprep(time =c(NA,rep("Tstop",3)), status =c(NA,paste0("stat",1:3)), data = LMall, 
                   keep = c(idOrig,"entry","LM",covs), trans = tmatCompete)
  attr(LMlong, "trans") <- tmatCompete
  class(LMlong) <- c("msdata", "data.table")
  #LMlong$Tstart<-LMlong$entry
  LMlong$Tstart<-LMlong$LM
  LMlong$time<-NULL
  LMlong <- expand.covs(LMlong, covs, longnames = F)
  return(LMlong)
}



pt.pred<-function(expo,pts,LM,pstats=2:4){
  pt<-pts[pts$exposure==expo,]
  ps<-matrix(as.numeric(cbind(pt[,paste0("pstate",pstats)])),ncol=1)
  data.frame(LM=rep(LM,3),
             exposure=rep(expo,3),
             trans=1:3,
             CIF=ps)
}

hr.beta<-function(x,cov_values){
  res=numeric(length(cov_values))
  for(i in 1:length(cov_values)){
    res[i] <-signif(x$coef[cov_values[i],"exp(coef)"], digits=3);#exp(beta)
  }
  return(res)
}

# calculate transition probabilities from msfits
trans.probs.LM<-function(LMlong,f, LMp){
  f<-as.formula(f)
  probs<-NULL
  ci.betas<-NULL
  for(landmark in LMp){
    c1 <- coxph(f, data = subset(LMlong,LM==landmark), method = "breslow")
    covs<-c("exposure1.1","exposure1.2","exposure1.3","exposure2.1",
            "exposure2.2","exposure2.3")
    HR.coefs<-hr.beta(x=summary(c1),cov_values=covs) 
    ci.betas<-rbind(ci.betas,
                    data.table(LM=rep(landmark,length(HR.coefs)),
                               exposure=covs,
                               HR=HR.coefs))
    
    msf.never <- msfit(c1, never, trans = tmatCompete)
    msf.cur <- msfit(c1, cur, trans = tmatCompete)
    msf.prev <- msfit(c1, prev, trans = tmatCompete)
    
    pts<-cbind(exposure=0:2,
               rbind(tail(probtrans(msf.never, predt=landmark,variance = F)[[1]],n=1),
                     tail(probtrans(msf.cur, predt=landmark,variance = F)[[1]],n=1),
                     tail(probtrans(msf.prev, predt=landmark,variance = F)[[1]],n=1)))
    
    data.table(do.call("rbind",lapply(0:2, pt.pred,pts=pts,LM=landmark)))
    probs<-rbind(probs,data.table(do.call("rbind",lapply(0:2, pt.pred,pts=pts,LM=landmark))))
  }
  return(list(ci.betas,probs))
}


## skew normal function for pacioli
rsn<-function (n = 1, xi = 0, omega = 1, alpha = 0, tau = 0, dp = NULL) 
{
  if (!is.null(dp)) {
    if (!missing(alpha)) 
      stop("You cannot set both 'dp' and component parameters")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    tau <- if (length(dp) > 3) 
      dp[4]
    else 0
  }
  if (tau == 0) {
    u1 <- rnorm(n)
    u2 <- rnorm(n)
    id <- (u2 > alpha * u1)
    u1[id] <- (-u1[id])
    z <- u1
  }
  else {
    delta <- alpha/sqrt(1 + alpha^2)
    truncN <- qnorm(runif(n, min = pnorm(-tau), max = 1))
    z <- delta * truncN + sqrt(1 - delta^2) * rnorm(n)
  }
  y <- as.vector(xi + omega * z)
  attr(y, "family") <- "SN"
  attr(y, "parameters") <- c(xi, omega, alpha, tau)
  return(y)
}


