##########################################################################################
# Plots and tables of the estimated hazard and cumulative incidences from the simulation ##
##########################################################################################

wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(wd,"/GlobalVars.R"))



library(reshape2)
library(scales) # for color names
library(RColorBrewer)
library(ggplot2)
library(xtable)


#################################################
## Set up graphical parameters and scenarios ####
#################################################
clrs<-hue_pal()(3)
clrs.paired<-brewer.pal(3, "Paired")
my_blue = brewer.pal(n = 6,"Blues")[c(5,1,6)] 
my_pastel = brewer.pal(n = 6,"Pastel1") 
clrs2=brewer.pal(n = 6,"BrBG")[c(2,1,3:6)]

theme_set(theme_minimal(base_size=14.5))



## labels for plots
state.lbls<-state.lbls2<-paste("State",1:6-1)
event.lbls<-state.lbls[4:6]
exposure.lbls<-state.lbls[1:3]
hr.lbls<-paste("Exposure",c(1:2)," vs. 0")
hr.labs<-hr.lbls
names(hr.labs)<-1:2
LT.lbls<-c("SN","Exp(0.1)","Exp(0.5)","Exp(1.5)")
LT.labs <- paste("LT",LT.lbls)
names(LT.labs) <- c(1:4)
exposure.labs <- hr.lbls
names(exposure.labs) <- 1:2
trans.labs<-event.lbls
names(trans.labs)<-1:3
interaction.lbls<-c(paste("1",c("constant","linear","qudadratic")), 
                    paste("2",c("constant","linear","qudadratic")))
names(interaction.lbls)<-c(1,1.1,1.2,2,2.1,2.2)
event.labs<-event.lbls
names(event.lbls)<-1:3


scenarios<-1:3
trans<-1:3


##################################
## Load the simulation results ###
##################################

Nparallel<-10 # number of parallel runs
nitPerRun<-1000
nit<-Nparallel*nitPerRun




simDataTranslist<-NULL
lmDataTranslist<-NULL
probsCrudelist<-NULL
coefsCrudelist<-NULL
coefsSupermodelAlllist<-NULL
coefsSupermodellist<-NULL
probsSupermodellist<-NULL
msprobslist<-NULL



## combine results from the parallel simulation runs
for(nparallel in 1:Nparallel){
  ## combine results from iterations in one run
  reslist<-readRDS(paste0(results.path,"reslist",nparallel,".Rds"))
  simDataTranslisti<-lapply(reslist,function(l)l$simDataTranslist)
  lmDataTranslisti<-lapply(reslist,function(l)l$lmDataTrans)
  probsCrudelisti<-lapply(reslist,function(l)l$probsCrude)
  coefsCrudelisti<-lapply(reslist,function(l)l$coefsCrude)
  coefsSupermodelAlllisti<-lapply(reslist,function(l)l$coefsSupermodelAll)
  coefsSupermodellisti<-lapply(reslist,function(l)l$coefsSupermodel)
  probsSupermodellisti<-lapply(reslist,function(l)l$probsSupermodel)
  msprobslisti<-lapply(reslist,function(l)l$msprobs)
  
  ## append the results to a list
  simDataTranslist<-c(simDataTranslist,simDataTranslisti)
  lmDataTranslist<-c(lmDataTranslist,lmDataTranslisti)
  probsCrudelist<-c(probsCrudelist,probsCrudelisti)
  coefsCrudelist<-c(coefsCrudelist,coefsCrudelisti)
  coefsSupermodelAlllist<-c(coefsSupermodelAlllist,coefsSupermodelAlllisti)
  coefsSupermodellist<-c(coefsSupermodellist,coefsSupermodellisti)
  probsSupermodellist<-c(probsSupermodellist,probsSupermodellisti)
  msprobslist<-c(msprobslist,msprobslisti)
}


lmDataTrans<-rbindlist(lmDataTranslist,idcol="it")
coefsCrude<-rbindlist(coefsCrudelist,idcol="it")

probsCrude<-rbindlist(probsCrudelist,idcol="it")
probsSuper<-rbindlist(probsSupermodellist,idcol="it")
probsEtm<-rbindlist(msprobslist,idcol="it")
probsTrue<-readRDS(paste0(data.path,"probsTrue.Rds"))



###############################################################################################
## Simulation data description ##############################################################

for(sc in 1:3){
  cat("Scenario",sc)
  simDataSclist<-lapply(simDataTranslist, function(l) l[[sc]])
  simDatasc<-Reduce("+",simDataSclist)
  sum(simDatasc[,3:5])/nit
  
  events.sim<-round(simDatasc/nit,digits = 0)
  colnames(events.sim)<-state.lbls2[2:6]
  rownames(events.sim)<-state.lbls2[1:3]
  print(events.sim)
  print(xtableFtable(ftable(events.sim),caption =paste0("Rounded average sampled number of transition events in scenario ",sc,"."),label =paste("tab:trans_events Simu",sc)))
}

###############################################################################################
## Landmark analysis description ##############################################################

head(lmDataTrans)
## number at risk at landmarks -----------------------------------------------------------
LMlongavg<-lmDataTrans[,list("value"=sum(value)/nit),by=c("scenario","LM","variable","type")]


sc.labs <- paste("Scenario",1:3)
names(sc.labs) <- c(1:3)
type.labs <- c("Size of risk set","Outomes at horizon")
names(type.labs) <- c("nrisk","nevent")

plotrisksets<-ggplot(LMlongavg[scenario%in%1:2], aes(x = LM, y = value, fill = variable)) + 
  geom_bar(stat="identity")+
  labs(y="Average frequency")+
  scale_fill_manual(values=clrs2,name="State", labels=state.lbls2)+
  theme(legend.position="bottom", legend.box = "horizontal")+
  ggtitle(paste("Composition of landmark data"))+
  theme(strip.background = element_rect(fill="lightgrey"))+
  facet_grid(scenario~type ,labeller = labeller(scenario=sc.labs,  type=type.labs))
plotrisksets
 

###############################################################################################
## Estimates description ##############################################################
## Separate landmark estimates #------------------------------------------------------------
coefsCrude$trans<-substr(coefsCrude$exposure, nchar(coefsCrude$exposure) - 1, nchar(coefsCrude$exposure))
coefsCrude$exposure<-substr(coefsCrude$exposure,  1, nchar(coefsCrude$exposure)-2)
coefsCrude$trans<-sub("\\.", "", coefsCrude$trans)
coefsCrude[,"exposure":=factor(coefsCrude$exposure, labels=1:2)]
coefsCrude.ci<-coefsCrude[,list(HR=mean(HR),lower=quantile(HR,probs=0.025,na.rm = T),
                                upper=quantile(HR,probs=0.975,na.rm = T)),by=list(LM,exposure,trans,scenario)]

pd2<-position_dodge(0.3)

hr.true<-data.table(HRTrue=c(1,1,5,10/3,2,1.5),
                    exposure=rep(c(1,2),3),
                    trans=rep(1:3,each=2))
hr.true[,c("exposure","trans"):=list(factor(exposure),as.character(trans))]
coefsCrude<-merge(coefsCrude.ci,hr.true,by=c("exposure","trans"), all.x = T)
head(coefsCrude)

##################################################################
## Figure S3 Estimated hazard rations in crude landmark models ###
##################################################################

plHRs<-ggplot(data=coefsCrude, aes(x = LM, y = HR, colour=as.factor(scenario))) +
  geom_hline(aes(yintercept=1),color="darkgrey", linetype="solid",linewidth=1)+
  geom_line(aes(x=LM,y=HRTrue),color="black", linetype="dashed",linewidth=1)+
  geom_errorbar(aes(ymax =upper, ymin = lower), width=0, position = pd2)+
  geom_point(size=1.8,position = pd2)+
  scale_color_manual(values=clrs.paired,name="Scenario", labels=scenarios)+
  theme(legend.position = c(0.1, 0.8))+
  xlab("Landmark")+
  ylab("Hazard ratio")+
  theme(legend.position="bottom", legend.box = "horizontal")+
  facet_grid(trans~exposure,
             labeller = labeller(trans=trans.labs,exposure=hr.labs))+
  scale_x_continuous(breaks=unique(coefsCrude$LM))+
  theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="lightgrey"), panel.border = element_rect(color = "black",fill=NA),
        legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal",
        text = element_text(size = 16)
  )
plHRs
ggsave(paste0(wdimages,"S3_HRCrudeSimulation.tiff"),plHRs,width=8,height=6,dpi=300)

## Supermodel #------------------------------------------------------------
coefs<-rbindlist(coefsSupermodellist,idcol="it")

coefs$trans<-substr(coefs$variable, nchar(coefs$variable) - 1, nchar(coefs$variable))
coefs$variable<-substr(coefs$variable,  1, nchar(coefs$variable)-2)
coefs$trans<-sub("\\.", "", coefs$trans)
coefs[,"HR":=`exp(coef)`]

coefsSuper.ci<-coefs[,list(coef=mean(HR),lowerSuper=quantile(HR,probs=0.025),
                           upperSuper=quantile(HR,probs=0.975)),
                     by=list(variable,trans,scenario)]
coefsSuper.ciExpo<-coefsSuper.ci[substr(variable,1,4)=="expo"]
setnames(coefsSuper.ciExpo,old="coef",new="HRSuper")
coefsSuper.ciExpo[,c("exposure"):=list(factor(substr(variable,nchar(variable),nchar(variable))))]
coefsSuper.ciExpo[,"variable":=NULL]
coefs.combined<-merge(coefsCrude,coefsSuper.ciExpo,by=c("exposure","trans","scenario"), all.x = T)





################################################################################
## Table 5 Estimated hazard rations and interactions in landmark supermodels ###
################################################################################

coefsAll<-rbindlist(coefsSupermodelAlllist,idcol="it")

coefsAll$trans<-substr(coefsAll$variable, nchar(coefsAll$variable) - 1, nchar(coefsAll$variable))
coefsAll$variable<-substr(coefsAll$variable,  1, nchar(coefsAll$variable)-2)
coefsAll$trans<-sub("\\.", "", coefsAll$trans)
head(coefsAll)
coefsAll[,"HR":=`exp(coef)`]
coefs.all.ci<-coefsAll[,list(coef=mean(HR),lower=quantile(HR,probs=0.025),upper=quantile(HR,probs=0.975)),by=list(variable,trans,scenario)]


coefslist<-lapply(scenarios,function(sc){
  resCox<-coefsSuper.ci[scenario==sc]
  #print(resCox)
  temp<-apply(round(resCox[,4:6],digits=3),1, function(i) paste0(i[1]," (",i[2],",",i[3],")"))
  resCox$variable
  tabprintadj<-(t(data.table("Expsoure 1"=temp[1:3],
                             "Expsoure 2"=temp[4:6],
                             "$\\gamma_1$"=temp[7:9],
                             "$\\gamma_2$"=temp[10:12])))
  out<-data.table("Parameter"=rownames(tabprintadj),tabprintadj)
  return(out)
})

coefstab<-rbindlist(coefslist,idcol="Scenario")
colnames(coefstab)<-c("Scenario","Parameter",event.lbls)
temp<-xtable(coefstab, 
             caption=paste("Estimated regression parameters of the landmark supermodel in the simulated data",
                           "on an exponential scale. $95\\%$ Confidence Interval in brackets."), 
             label=paste("tab:supermodel coefficients"))
print(temp,sanitize.text.function=function(x){x},include.rownames = F)









#############################################################
## Figure 4, Figure S4 Estimated transition probabilities ###
#############################################################
probs.combined<-rbindlist(list("Crude LM"=probsCrude,"Super LM"=probsSuper, "Mstate"=probsEtm),idcol="Model")
probs.combined[,"Exposure":=factor(paste("Exposure",exposure))]
probs.combined[,"Scenario":=factor(scenario,levels=1:3,labels = paste("Scenario",1:3))]
probs.combined[,"Model":=factor(Model,levels=c("Mstate","Crude LM", "Super LM"))]
probsTrue[,"Exposure":=factor(paste("Exposure",exposure))]
probs.combined.estimates<- probs.combined[,list(
  CIF=mean(CIF,na.rm=T),lower=quantile(CIF,probs=0.025,na.rm=T),upper=quantile(CIF,probs=0.975,na.rm=T)), 
  by=list(Model,LM,Exposure,trans,Scenario)]



pd2<-position_dodge(0.9)
legend.name="Model"

# Figure 4 (SAB)
plProbsSAB<-ggplot(probs.combined.estimates[trans==2&LM%in%seq(1,13,by=2)], aes(x = LM, y = CIF, colour=Model)) +
  geom_line(data=probsTrue[trans==2],aes(x=LM,y=CIF),color="black", linetype="dashed",linewidth=1,alpha=0.5)+
  geom_errorbar(aes(ymax =upper, ymin = lower), width=0, position = pd2) +
  geom_point(position = pd2,size=2)+
  scale_x_continuous(breaks=seq(1,13,by=2),labels=seq(1,13,by=2))+
  ylab("Predicted cumulative incidence for state 4 at horizon")+
  xlab("Landmark")+
  facet_grid(as.factor(Scenario)~Exposure,
             labeller = labeller(scenario=trans.labs,exposure=exposure.labs))+
  theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="lightgrey"), panel.border = element_rect(color = "black",fill=NA),
        legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal",
        text = element_text(size = 16)
  )
plProbsSAB
ggsave(paste0(wdimages,"Unseld_Figure4.tiff"),plProbsSAB,width=8,height=8, unit="in",dpi=300)



# Figure S4 (ETOP)
plProbsETOP<-ggplot(probs.combined.estimates[trans==3&LM%in%seq(1,13,by=2)], aes(x = LM, y = CIF, colour=Model)) +
  geom_line(data=probsTrue[trans==3],aes(x=LM,y=CIF),color="black", linetype="dashed",linewidth=1,alpha=0.5)+
  geom_errorbar(aes(ymax =upper, ymin = lower), width=0, position = pd2) +
  geom_point(position = pd2,size=2)+
  scale_x_continuous(breaks=seq(1,13,by=2),labels=seq(1,13,by=2))+
  ylab("Predicted cumulative incidence for state 5 at horizon")+
  xlab("Landmark")+
  facet_grid(as.factor(Scenario)~Exposure,
             labeller = labeller(scenario=trans.labs,exposure=exposure.labs))+
  theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="lightgrey"), panel.border = element_rect(color = "black",fill=NA),
        legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal",
        text = element_text(size = 16)
  )
plProbsETOP
ggsave(paste0(wdimages,"S4_CIFsimulation_ETOP.tiff"),plProbsETOP,width=8,height=8, unit="in",dpi=300)





















