#############################################################
# Calculation of the true transition probabilities ##########
#############################################################

## Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(wd,"/GlobalVars.R"))
runfiles.path<-paste0(wd,"/MultipleScenarios/RunFiles/")
dataMultiple.path<-paste0(wd,"/MultipleScenarios/Data/")

wdimages<-paste(wdimages,"MultipleScenarios/") # save plots

source(paste0(runfiles.path,"AnaLandmarkScs.R"),local = F)
source(paste0(runfiles.path,"LMSimulationScs.R"),local = F)

library(gridExtra)

tmat<-transMat(x = list(c(1, 3:5)+1, c(2,3:5)+1,c(1,3:5)+1, c(), c(), c()))
Settings<-readRDS(paste0(runfiles.path,"Settings.Rds"))
Settings<-Settings[order(setting)]
ParsAllSettings<- Settings[!duplicated(Settings$setting),c("setting","HRSABCur","HRSABPrev",
                                                           "HRETOPCur","HRETOPPrev",
                                                           "RateCurPrev")]
ParsAllSettings

direction<-"forward"
variance<-F
covariance<-F
method<-"greenwood"


## Function for compuation of transition probabilities adapted from the 
# probtrans() function in the mstate package
Myprobtrans<-function (predt,Haz) {
  object<-list(trans=tmat,Haz=Haz)

  trans <- object$trans
  transit <- to.trans2(trans)
  numtrans <- nrow(transit)
  stackhaz <- object$Haz
  
  for (i in 1:numtrans) stackhaz$dhaz[stackhaz$trans == i] <- diff(c(0, 
                                                                     stackhaz$Haz[stackhaz$trans == i]))
  if (direction == "forward") stackhaz <- stackhaz[stackhaz$time > predt, ]

  untimes <- sort(unique(stackhaz$time))
  TT <- length(untimes)
  S <- nrow(trans)
  
  if (direction == "forward") {
    if (variance == TRUE) 
      res <- array(0, c(TT + 1, 2 * S + 1, S))
    else res <- array(0, c(TT + 1, S + 1, S))
    res[1, 1, ] <- predt
    for (j in 1:S) res[1, 1 + j, ] <- rep(c(0, 1, 0), c(j - 
                                                          1, 1, S - j))
    if (variance) 
      res[1, (S + 2):(2 * S + 1), ] <- 0
  }
  
  P <- diag(S)
 
  for (i in 1:TT) {
    idx <- ifelse(direction == "forward", i, TT + 1 - i)
    tt <- untimes[idx]
    Haztt <- stackhaz[stackhaz$time == tt, ]
    lHaztt <- nrow(Haztt)
    IplusdA <- diag(S)
    for (j in 1:lHaztt) {
      from <- transit$from[transit$transno == Haztt$trans[j]]
      to <- transit$to[transit$transno == Haztt$trans[j]]
      IplusdA[from, to] <- Haztt$dhaz[j]
      IplusdA[from, from] <- IplusdA[from, from] - Haztt$dhaz[j]
    }
    if (any(diag(IplusdA) < 0)) 
      warning("Warning! Negative diagonal elements of (I+dA); the estimate may not be meaningful. \n")

 if (method == "greenwood") {
      if (direction == "forward") {
        
        P <- P %*% IplusdA
      }
      
    }
   
    if (direction == "forward") {
      res[idx + 1, 1, ] <- tt
      res[idx + 1, 2:(S + 1), ] <- t(P)
      
    }
   
  }
  
  res2 <- vector("list", S)
  for (s in 1:S) {
    tmp <- as.data.frame(res[, , s])
    if (min(dim(tmp)) == 1) 
      tmp <- res[, , s]
    if (variance) 
      names(tmp) <- c("time", paste("pstate", 1:S, sep = ""), 
                      paste("se", 1:S, sep = ""))
    else names(tmp) <- c("time", paste("pstate", 1:S, sep = ""))
    res2[[s]] <- tmp
  }

  res2$trans <- trans
  res2$method <- method
  res2$predt <- predt
  res2$direction <- direction
  class(res2) <- "probtrans"
  return(res2)
}


# Help functions -------------------------------------------------------------
to.trans2<-function (trans) 
{
  dm <- dim(trans)
  if (dm[1] != dm[2]) 
    stop("transition matrix should be square")
  S <- dm[1]
  mx <- max(trans, na.rm = TRUE)
  res <- matrix(NA, mx, 3)
  res[, 1] <- 1:mx
  transvec <- as.vector(trans)
  for (i in 1:mx) {
    idx <- which(transvec == i)
    res[i, 2:3] <- c((idx - 1)%%S + 1, (idx - 1)%/%S + 1)
  }
  res <- data.frame(res)
  names(res) <- c("transno", "from", "to")
  res$from[res$from == 0] <- S
  statesfrom <- dimnames(trans)[[1]]
  if (is.null(statesfrom)) 
    statesfrom <- 1:S
  statesto <- dimnames(trans)[[2]]
  if (is.null(statesto)) 
    statesto <- 1:S
  res$fromname <- statesfrom[res$from]
  res$toname <- statesto[res$to]
  res$transname <- paste(res$fromname, res$toname, sep = " -> ")
  return(res)
}


##############################################################################################
# calculate transition probabilities for each setting from true cumulative hazards
probsTrueSettings<-function(ParsSetting){
  cum.haz<-specifyCumHaz(ParsSetting)
  
  pt0list<-lapply(LMps,Myprobtrans, Haz=cum.haz)
  ptslist<-lapply(pt0list,function(pt) {
    cbind(exposure=0:2,
          rbind(tail(pt[[1]],n=1), # from never exposed before
                tail(pt[[2]],n=1), # from currently exposed
                tail(pt[[3]],n=1)))}) # from previouxly exposed
  
  probslist<-lapply(1:length(LMps), function(i) {
    data.table(do.call("rbind",lapply(0:2, pt.pred,pts=ptslist[[i]],LM=LMps[i],pstats=4:6)))
  
  })
  
  probs<-do.call("rbind",probslist)

  return(probs)
}

##############################################################################################
LMps<-seq(0,13, by=0.05) # landmark points with distance as in cumulative hazard

#################################################################################################
### Transition probabilities for the first simulation -----------------------------------------
ParsSetting=head(ParsAllSettings[HRSABCur==5& HRSABPrev==10/3 &HRETOPCur==2  & RateCurPrev==1,],n=1)
ParsSetting[,"HRETOPPrev"]<-1.5

## cacluate probabilities
probsTrue<-probsTrueSettings(ParsSetting)
saveRDS(probsTrue,paste0(wd,"/Data/probsTrue.Rds"))





##############################################################################################
### Transition probabilities for the additional settings --------------------------------------

# calculate transition probabilities for each seting from true cumulative hazards
probsTrueSettings<-function(setti){
  ParsSetting<-ParsAllSettings[setting==setti]
  # print(ParsSetting)
  cum.haz<-specifyCumHaz(ParsSetting)
  
  pt0list<-lapply(LMps,Myprobtrans, Haz=cum.haz)
  ptslist<-lapply(pt0list,function(pt) {
    cbind(exposure=0:2,
          rbind(tail(pt[[1]],n=1), # from never exposed before
                tail(pt[[2]],n=1), # from currently exposed
                tail(pt[[3]],n=1)))}) # from previouxly exposed
  
  probslist<-lapply(1:length(LMps), function(i) {
    data.table(do.call("rbind",lapply(0:2, pt.pred,pts=ptslist[[i]],LM=LMps[i],pstats=4:6)))
    
  })
  
  probs<-do.call("rbind",probslist)
  
  return(probs)
}


## cacluate probabilities
settings<-levels(Settings$setting)
probsTruelist<-lapply(settings,probsTrueSettings)
probsTrue<-rbindlist(probsTruelist,idcol="setting")
saveRDS(probsTrue,paste0(dataMultiple.path,"probsTrueScs.Rds"))













