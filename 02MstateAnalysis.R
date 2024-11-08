## Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)

source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"/HelpFunctionsAnaGyrase.R"), local = T)

## Load packages and functions
library(mvna)


dMSM<-readRDS(paste0(data.path,"dMSM_anonymized.Rds"))
dMSM2<-readRDS(paste0(data.path,"dMSM2.Rds"))
dMSMexpand<-readRDS(paste0(data.path,"dMSMexpand.Rds"))



#############################################
## Test for independent left-truncation #####
#############################################
reslist<-lapply(3:5, function(i) data.table(confint.beta(summary(coxph(Surv(entry, exit, to==i) ~ study.entry+from+
                                                                         cluster(idChild), dMSM2)),
                                                         cov_values = "study.entry")))
# adjustment for age group
reslistadj<-lapply(3:5, function(i) data.table(confint.beta(summary(coxph(Surv(entry, exit, to==i) ~ study.entry+from+age.gr+
                                                                            cluster(idChild), dMSM2)),
                                                            cov_values = "study.entry")))
# all-cause event
resAll<-rbind(data.table(confint.beta(summary(coxph(Surv(entry, exit, to%in%3:6) ~ study.entry+from+cluster(idChild), dMSM2)),
                                      cov_values = "study.entry")),
              data.table(confint.beta(summary(coxph(Surv(entry, exit, to%in%3:6) ~ study.entry+from+age.gr+cluster(idChild), 
                                                    dMSM2)),cov_values = "study.entry")))

templ<-list(rbindlist(reslist,idcol="outcome"),rbindlist(reslistadj,idcol="outcome"))
leftTruncTests<-rbindlist(templ,idcol="model")
leftTruncTests[,c("model","outcome"):=list(ifelse(model==1,"Unadjusted","Adjusted"),
                                           ifelse(outcome==1,"Life birth",ifelse(outcome==2,"SAB","ETOP")))]
leftTruncTests<-rbind(leftTruncTests,cbind(model=c("Unadjusted","Adjusted"),outcome="All-cause",resAll))
leftTruncTests<-leftTruncTests[order(model, decreasing = T)]
leftTruncTests
print(xtable(leftTruncTests,caption=paste("Coefficient estimates of the left-truncation effect in Cox models for the pregnancy outcomes.",
                                    "All-cause, time to any pregnancy outcome."), label = "tab: independent left-truncation Test", digits=4),
      include.rownames = F)




############################
## Non-parametric model ###
############################
mvnaConfint<-function(mvnaObj, tr.choice){
  object <- mvna::summary.mvna(mvnaObj, level = 0.95, var.type ="aalen", 
                               ci.fun ="log")[tr.choice]
  naeslist<-lapply(object,function(x) data.table("time"=x$time, "NAE"=x$na,"lower"=x$lower,"upper"=x$upper))
  naes<-rbindlist(naeslist,idcol="trans")
  return(naes)
}

## NAE in multi-state model 1 with currently and previously exposed state -----------------
state.lbls<-0:6
tmat<-!is.na(transMat(x = list(c(1, 3:6)+1, c(2,3:6)+1,c(1,3:6)+1, c(), c(), c(),c()), 
                      names = state.lbls))
tmat
setnames(dMSM,old="idChild",new="id")
mvna.gyrase<-mvna(dMSM,state.lbls,tmat,cens=NULL)

## NAE in multi-state model 2 with exposed before state -----------------------------------
state.lbls2<-c(0,1,3:6)
tra
setnames(dMSM2,old="idChild",new="id")
mvna.gyrase2<-mvna(dMSM2,state.lbls2,tra,cens=NULL)




# table with number of events and risk set at a selective samlple of times
dt<-data.table(time=mvna.gyrase2$time,
               "0_1"=cumsum(mvna.gyrase2$n.event[1,2,]),
               "0_3"=cumsum(mvna.gyrase2$n.event[1,3,]),
               "1_3"=cumsum(mvna.gyrase2$n.event[2,3,]),
               "0_4"=cumsum(mvna.gyrase2$n.event[1,4,]),
               "1_4"=cumsum(mvna.gyrase2$n.event[2,4,]),
               "0_5"=cumsum(mvna.gyrase2$n.event[1,5,]),
               "1_5"=cumsum(mvna.gyrase2$n.event[2,5,]),
               "0_6"=cumsum(mvna.gyrase2$n.event[1,6,]),
               "1_6"=cumsum(mvna.gyrase2$n.event[2,6,]),
               "0"=mvna.gyrase2$n.risk[,1],
               "1"=mvna.gyrase2$n.risk[,2])

dtLM<-dt[time%in%c(1:13,15,20,23,25,30,35,40,42,max(dt$time))]
temp<-diff(as.matrix(dtLM[,2:10]))
dtLMdiff<-data.table(time=dtLM$time,rbind(rep(0,ncol(temp)),temp),dtLM[,11:12])
dtLMdiff


## transition probabilities for outcomes at horizon ###############
mspred<-function(tti,msf0){
  pt0 <- probtrans(msf0, direction="forward", predt = tti, method = "greenwood",variance = T)
  pts<-cbind(exposure=0:1,
             rbind(tail(pt0[[1]],n=1), # from never exposed before
                   tail(pt0[[2]],n=1))) # from exposed before
  probs<-data.table(do.call("rbind",lapply(0:1, pt.pred,pts=pts,LM=tti, pstats=3:5)))
  return(probs)
} 

LMps<-5:13
c0<-coxph(Surv(entry, exit, status)~strata(trans)+cluster(idChild),data=dMSMexpand)
msf0 <- msfit(object = c0, vartype = "greenwood", trans =tra.num)
probslist<-lapply(LMps, mspred, msf=msf0)
probs<-rbindlist(probslist)
probs
saveRDS(probs,file=paste0(wdresults,"probsMstate.Rds"))



















