## Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)

source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"/HelpFunctionsAnaGyrase.R"), local = T)



is<-as.character(3:5) # Competing outcome events of interest


dMSM2<-readRDS(paste0(data.path,"dMSM2.Rds"))




########################################################################################
## Exposure as a time-dependent covariate ##############################################
########################################################################################

## Prepare covariates
dMSM2[,"from":=factor(from)]
dMSM2[,"age.gr":=factor(age.gr,levels=c("20-35 years","<=20 years",">35 years"))]


# Cox models for each outcome event
coxTDC<-function(cov_values="from"){
  foreach(i=iter(is),.combine = rbind)%do%{
    # robust confidence intervals due to multiple births and multiple transitions of some women
    formula_eventi<-as.formula(paste0("Surv(entry, exit, to == ",i,")~ ",paste(cov_values,collapse="+"),"+cluster(idWoman)")) 
    tab_tdc_eventi<-summary(coxph(formula_eventi, data = dMSM2, method = "breslow"))
    out<-data.table(confint.beta(tab_tdc_eventi,"cov_values"=rownames(tab_tdc_eventi$coefficients)))
    out[,"Variable":=rownames(tab_tdc_eventi$coefficients)]
    out[,"Outcome":=i]
    return(out)
  }
}


## Table 1: Cox Model with exposure as time-dependent covariate ###
# (estimates don't exactly match those from the manuscript because of the anonymization)

rescoxTDCundadj<-coxTDC(cov_values=c("from"))
rescoxTDCadj<-coxTDC(cov_values=c("from","age.gr"))
rescoxTDC<-rbindlist(list(rescoxTDCundadj,rescoxTDCadj),idcol="Adjustment",fill=TRUE)

rescoxTDC[,"Outcome":=factor(Outcome,labels=event.lbls[-4])]
rescoxTDC[,"Adjustment":=factor(Adjustment,labels=c("Unadjusted","Adjusted"))]
rescoxTDC
saveRDS(rescoxTDC,file=paste0(wdresults,"GyrasecoxTDC.Rds"))


tabcoxTDC<-copy(rescoxTDC)
tabcoxTDC[,"HR":=apply(rescoxTDC, 1,function(tabi) paste0(tabi["HR"]," (",tabi["lower .95"],",",tabi["upper .95"],")"))]
tabcoxTDC[,"Variable":=factor(Variable,
                                  levels=c("from1","age.gr<=20 years","age.gr>35 years"),
                                  labels = c("Fluoroquionolone", "Age $\\leq 20$ years", "Age $> 35$ years"))]
setorder(tabcoxTDC,-Variable)


tabcoxTDC
tabcoxTDCwide<-data.table::dcast(tabcoxTDC,Variable~Outcome+Adjustment,value.var="HR",sep = "__")
tabcoxTDCwideadj<-data.table::dcast(tabcoxTDC[Adjustment=="Adjusted"],Variable~Outcome,value.var="HR")
print(xtable(tabcoxTDCwideadj,  caption = paste("Estimated hazard ratios and 95\\% pointwise robust confidence intervals for experiencing the different pregnancy outcomes (live birth, spontaneous abortion (SAB), elective",
                                             "termination of pregnancy (ETOP)in Cox models incorporating fluoroquinolone exposure (exposed vs. never exposed before) as time-dependent covariate",
                                             "and age at study entry as baseline covariate."), 
                      label="tab:gyrase_tdc"), include.rownames = F, sanitize.text.function = function(x)x)



