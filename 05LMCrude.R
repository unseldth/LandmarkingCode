## Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"HelpFunctionsAnaGyrase.R"), local = T)


LMlong<-readRDS(paste0(data.path,"/LMlong.Rds"))

############################################################################################
## Crude landmark model ####################################################################
############################################################################################
LMlong<-LMlong%>%filter(to!=6) # exclude stillbirths because there are too few events


# unadjusted ################
nexposure<-data.frame(exposure = factor(rep(0, nevents),levels=0:1), trans=1:nevents)
nexposure<-expand.exposure(nexposure, "exposure")

before.expo<-data.frame(exposure = factor(rep(1, nevents),levels=0:1), trans=1:nevents)
before.expo<-expand.exposure(before.expo, "exposure")

covsexposure<-paste0("exposure.",1:(nevents-1))


f<-paste0('Surv(Tstart,Tstop, status) ~',paste(covsexposure,collapse="+"),'+strata(trans)+cluster(idWoman)')
res.ever<-trans.probs.LM2(f) 
coefsCrude<-res.ever[[1]]
probsCrude<-res.ever[[2]]
saveRDS(probsCrude,file = paste0(wdresults,"probsCrude.Rds"))

# adjusted ################
# data frame with reference categories and variable exposure level
expo.base.covs<-function(expo) {
  data.frame(exposure = factor(rep(expo, nevents),levels=0:1),
             age.gr=factor(rep("20-35 years", nevents),levels=age.lbls),
             trans=1:nevents)
}

nexposure<-expo.base.covs(expo=0)
nexposure<-expand.exposure(nexposure, c("exposure","age.gr"))
covsadj<-c(paste0("age.gr1.",1:(nevents-1)),paste0("age.gr2.",1:(nevents-1)))

before.expo<-expo.base.covs(expo=1)
before.expo<-expand.exposure(before.expo, c("exposure","age.gr"))


f.adj<-paste0('Surv(Tstart,Tstop, status) ~',paste(c(covsexposure,covsadj),collapse="+"),'+strata(trans)+cluster(idWoman)')


res.adj<-trans.probs.LM2(f.adj,cov_values = c(covsexposure,covsadj)) 
coefsCrude.adj<-res.adj[[1]]
probsCrude.adj<-res.adj[[2]]



#####################################################################
## Table 2 Hazard Ratios for the crude landmark model ###############
#####################################################################

coefslist<-list(coefsCrude,coefsCrude.adj)
coefs.ests<-rbindlist(coefslist,idcol="Model")
coefs.ests[,"Model":=models.lbls[Model]]
coefs.ests$trans<-substr(coefs.ests$Variable, nchar(coefs.ests$Variable) - 1, nchar(coefs.ests$Variable))
coefs.ests$Variable<-substr(coefs.ests$Variable,  1, nchar(coefs.ests$Variable)-2)
coefs.ests$trans<-sub("\\.", "", coefs.ests$trans)
setnames(coefs.ests, old=c("lower .95","upper .95" ), new=c("lower","upper"))
coefs.ests[,"Variable":=factor(Variable,levels=c("exposure","age.gr1","age.gr2"),
                               labels=c("Fluoroquinolone","Age $\\leq$20 years", "Age $>$35 years"))]
coefs.ests[,"Outcome":=factor(trans,labels=event.lbls[-4])]
setnames(coefs.ests,old="Model",new="Adjustment")
coefs.ests

temp1<-as.data.table(copy(coefs.ests[trans!=4]))
temp1[,"Model":="Crude LM"]




coefs.ests[lower==0]
coefs.ests2<-coefs.ests%>%dplyr::mutate(
  "Landmark"=LM,
  "Adjustment"=factor(Adjustment,levels=c("Unadjusted","Adjusted"),labels=c("Unadjusted","Age-adjusted")),
  "HR (95\\% CI)"=ifelse(lower==0,"-",sprintf("%.3f (%.3f,%.3f)", HR,lower,upper))
)
coefswide<-dcast(coefs.ests2,Variable+Landmark~Outcome+Adjustment,value.var ="HR (95\\% CI)" )
coefswide

coefswideadj<-dcast(coefs.ests2[Adjustment=="Age-adjusted"],Variable+Landmark~Outcome,value.var ="HR (95\\% CI)" )
coefswideadj

xtab_coefsadj<-xtable(coefswideadj,
                      caption=paste("Estimated hazard ratios (HR) and 95\\% pointwise robust confidence intervals in",
                                    "crude landmark models for live birth, spontaneous abortion (SAB), and",
                      "elective termination of pregnancy (ETOP). Fluoroquinolone: HR of prior exposure vs. no prior",
                      "exposure. Age: HRs of the lower age group ($\\leq$ 20 years) or the older age group ($>$ 35 years)",
                      "vs. the middle age group (20$-$35 years)."),
                      label = "tab:hrgyrasecrude")
print(xtab_coefsadj,sanitize.text.function = function(x) x, include.rownames = F)





