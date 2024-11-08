#################################################################################################################
# Plots and tables of the estimated hazard and cumulative incidences from the simulation in additional settings ##
#################################################################################################################
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
runfiles.path<-paste0(wd,"/RunFiles/") # to load cumHazSimu.Rds and for LMSimulation.R (which calls AnaLandmark.R)
data.path<-paste0(wd,"/Data/")  


probsTrue<-readRDS(paste0(data.path,"probsTrueScs.Rds"))
Settings<-readRDS(paste0(runfiles.path,"Settings.Rds"))


library(reshape2)
library(scales) 
library(RColorBrewer)
library(ggplot2)
library(xtable)
library(gridExtra)

source(paste0(wd,"/../GlobalVars.R"))

wdimages<-paste0(wd,"/../../Figures/Simulation/MultipleScenarios/")  

#################################################
## Set up graphical parameters and scenarios ####
#################################################
n<-6000 # number of women per iteration


settings<-1:nrow(Settings)
nsettings<-nrow(Settings) ## change
nLTs<-4
trans<-1:3

clrs<-hue_pal()(nsettings)
clrs.paired<-brewer.pal(nsettings, "Paired")
my_blue = brewer.pal(n = 3,"Blues")[c(2,3,1)] 

## labels for plots
LT.lbls<-c("SN","Exp(0.1)","Exp(0.5)","Exp(1.5)")
LT.labs <- paste("LT",LT.lbls)
names(LT.labs) <- c(1:4)
exposure.labs <- hr.lbls
names(exposure.labs) <- 1:2
setting.labs <- paste("Setting",1:nsettings)
names(setting.labs) <- c(1:nsettings)
trans.labs<-event.lbls
names(trans.labs)<-1:3
interaction.lbls<-c(paste("Currently",c("constant","linear","qudadratic")), 
                    paste("Previously",c("constant","linear","qudadratic")))
names(interaction.lbls)<-c(1,1.1,1.2,2,2.1,2.2)
event.labs<-event.lbls
names(event.lbls)<-1:3

## append true probabilities for scenario 9
probsTrue9<-probsTrue[setting==6]
probsTrue9[,"setting":=9]
probsTrue<-rbind(probsTrue,probsTrue9)
probsTrue[,"setting":=factor(setting,levels=1:nsettings)]
print(xtable(Settings[,c('setting','HRSABCur','HRSABPrev','HRETOPCur','HRETOPPrev','RateCurPrev','pStartNever')], 
             caption="Hazard modifications for the additional settings",
             label="tab: Multiple scenarios hazards modifications"), include.rownames=F)

##################################
## Load the simulation results ###
##################################
simDataTrans<-readRDS(paste0(results.path,"simDataTrans.Rds"))
lmDataTrans<-readRDS(paste0(results.path,"lmDataTrans.Rds"))
msprobs<-readRDS(paste0(results.path,"msprobs.Rds"))
coefsTdc<-readRDS(paste0(results.path,"coefsTdc.Rds"))
coefsCrude<-readRDS(paste0(results.path,"coefsCrude.Rds"))
probsCrude<-readRDS(paste0(results.path,"probsCrude.Rds"))
coefsSupermodelAll<-readRDS(paste0(results.path,"coefsSupermodelAll.Rds"))
coefsSupermodel<-readRDS(paste0(results.path,"coefsSupermodel.Rds"))
probsSupermodel<-readRDS(paste0(results.path,"probsSupermodel.Rds"))


head(simDataTrans)
head(lmDataTrans)
head(msprobs)
head(coefsTdc)
head(coefsCrude)
  
###############################################################################################
## Simulation data description ##############################################################
## Number of transitions ###
setnames(simDataTrans,old=as.character(2:6),new=paste0("To",2:6))
nit<-length(unique(simDataTrans$iteration))
events.sim<-simDataTrans[,list('Currently'=sum(To2)/nit,
                             'Previously'=sum(To3)/nit,
                             'Life Birth'=sum(To4)/nit,
                             'SAB'=sum(To5)/nit,
                             'ETOP'=sum(To6)/nit),
                             by=c("setting", "LT","from")]
## Table
x<-events.sim[order(from,LT)]
x[,c("from","setting"):=list(exposure.lbls[as.numeric(from)],factor(setting))]
print(xtable(x,caption="Rounded average number of transitions over 1000 iterations.",
             label = "tab:nTrans multiple scenarios", digits=0), include.rownames = F,
      size="\\fontsize{10pt}{10pt}\\selectfont", tabular.environment = "longtable", floating = F)


## Plots
events.sim.m<-as.data.table(melt(events.sim,id=c("setting", "LT","from")))
events.sim.m[,c("setting"):=list(factor(setting))]

LTsLayout=matrix(1:4,nrow=2,byrow=T)
ymaxmat<-events.sim.m[,list("ymaxfrom"=max(value)),by=c("from")]
plotEvents<-function(LTis,fromi){
  if(LTis==1){
    pl1<-ggplot(events.sim.m[from==fromi & LT %in% LTsLayout[LTis,]],aes(fill=setting,y=value,x=variable))+
      geom_bar(stat = 'identity', position = 'dodge')+
      facet_grid(.~LT, labeller = labeller(LT=LT.labs))+
      scale_y_sqrt(limits=c(0,ymaxmat[from==fromi]$ymaxfrom+50))+
      scale_fill_brewer(palette="Paired")+
      theme_minimal()+
      xlab("To")+
      ylab("Frequency (sqrt scale)")+
      theme_minimal()+
      theme(legend.position = "none")+
      theme(strip.background = element_rect(fill="lightgrey"))
    pl1
    return(pl1)
  }
  pl2<-ggplot(events.sim.m[from==fromi & LT %in% LTsLayout[LTis,]],aes(fill=setting,y=value,x=variable))+
    geom_bar(stat = 'identity', position = 'dodge')+
    facet_grid(.~LT, labeller = labeller(LT=LT.labs))+
    scale_y_sqrt(limits=c(0,ymaxmat[from==fromi]$ymaxfrom+50))+
    scale_fill_brewer(palette="Paired")+
    theme_minimal()+
    xlab("To")+
    ylab("Frequency (sqrt scale)")+
    theme_minimal()+
    theme(strip.background = element_rect(fill="lightgrey"))+
    theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")+
    guides(fill = guide_legend(nrow = 1))
  return(pl2)
}

for(fromi in 1:3){
  nevents.sim.pl<-lapply(1:2,plotEvents,fromi=fromi)
  ml <- grid.arrange(nevents.sim.pl[[1]],nevents.sim.pl[[2]], heights=c(10,12),
                    ncol=1,nrow=2,top=paste("Average transition numbers from",
                                                                         exposure.lbls[fromi],"exposed"))
  print(ml)
}



###########################################################
## Figure S5, S6: Landmark datasets description ###########
###########################################################
LMlong<-lmDataTrans
LMlong[,c("LT","count"):=list(factor(LT),value)]
nrisk<-LMlong[type=="nrisk"]
nevent<-LMlong[type=="nevent"]

nrisk.m<-nrisk[,list("count"=sum(count)/nit),by=c("setting","LT","LM","variable")]
nevent.m<-nevent[,list("count"=sum(count)/nit),by=c("setting","LT","LM","variable")]

plotrisksets<-function(sets){
  nrisk.pl<-ggplot(nrisk.m[setting%in%sets], aes(fill=variable, y=count, x=LM)) + 
    geom_bar(stat = 'identity', position = 'stack') +
    facet_grid(setting~ LT,labeller = labeller(LT=LT.labs,setting=setting.labs))+
    ggtitle(paste("Composition of landmark data"))+
    scale_fill_manual(values=my_blue,name="Exposure", labels=exposure.lbls)+
    xlab("Landmark time points (gestational week)")+
    labs(y="Average frequency")+
    theme_minimal()+
    theme(legend.position="bottom", legend.box = "horizontal")+
    scale_y_continuous(breaks=seq(0,6000,by=2000))+
    theme(strip.background = element_rect(fill="lightgrey"))+
    theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")+
    guides(fill = guide_legend(nrow = 1))
    print(nrisk.pl)
  ggsave(paste0(wdimages,"S5_nrisk.tiff"), nrisk.pl,width=8,height=10,dpi=300)
} 

plotrisksets(1:nsettings)

plotevents<-function(sets){
  nevent.pl<-ggplot(nevent.m[setting%in%sets], aes(fill=variable, y=count, x=LM)) + 
    geom_bar(stat = 'identity', position = 'stack') +
    facet_grid(setting~ LT,labeller = labeller(LT=LT.labs,setting=setting.labs))+
    scale_fill_brewer(palette="Paired",labels=event.labs)+
    xlab("Landmark time points (gestational week)")+
    ylab("Average frequency")+
    theme_minimal()+
    ggtitle("Outcomes at horizon")+
    theme(strip.background = element_rect(fill="lightgrey"))+
    theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")+
    guides(fill = guide_legend(nrow = 1,title = "Outcome"))
  return(nevent.pl)
}
 

plotevents(1:nsettings)
ggsave(paste0(wdimages,"S6_nevent.tiff"),width=8,height=10,dpi=300)


########################################################
## Figures S7-S24: coefficient and CIF estimates #######
########################################################
## TDC estimates #------------------------------------------------------------
coefsTdc.ci<-coefsTdc[,list(HR=mean(HR),lower=quantile(HR,probs=0.025,na.rm = T),upper=quantile(HR,probs=0.975,na.rm = T)),
                      by=list(exposure,trans,setting,LT)]
coefsTdc.ci[,c("exposure","trans"):=list(factor(exposure),as.character(trans))]


## Separate landmark estimates #------------------------------------------------------------
coefsCrude$trans<-substr(coefsCrude$exposure, nchar(coefsCrude$exposure) - 1, nchar(coefsCrude$exposure))
coefsCrude$exposure<-substr(coefsCrude$exposure,  1, nchar(coefsCrude$exposure)-2)
coefsCrude$trans<-sub("\\.", "", coefsCrude$trans)
coefsCrude[,"exposure":=factor(coefsCrude$exposure, labels=1:2)]
coefsCrude.ci<-coefsCrude[,list(HR=mean(HR),lower=quantile(HR,probs=0.025,na.rm = T),upper=quantile(HR,probs=0.975,na.rm = T)),
                          by=list(LM,exposure,trans,setting,LT)]

## Supermodel #------------------------------------------------------------
coefs<-coefsSupermodel

coefs$trans<-substr(coefs$variable, nchar(coefs$variable) - 1, nchar(coefs$variable))
coefs$variable<-substr(coefs$variable,  1, nchar(coefs$variable)-2)
coefs$trans<-sub("\\.", "", coefs$trans)

coefsSuper.ci<-coefs[,list(coef=mean(HR),lowerSuper=quantile(HR,probs=0.025),upperSuper=quantile(HR,probs=0.975)),by=list(variable,trans,setting,LT)]
coefsSuper.ciExpo<-coefsSuper.ci[substr(variable,1,4)=="expo"]
setnames(coefsSuper.ciExpo,old="coef",new="HRSuper")
coefsSuper.ciExpo[,c("exposure"):=list(factor(substr(variable,nchar(variable),nchar(variable))))]
coefsSuper.ciExpo[,"variable":=NULL]

coefs.combined<-merge(coefsCrude.ci,coefsTdc.ci,by=c("exposure","trans","setting","LT"), all.x = T)
coefs.combined<-merge(coefs.combined,coefsSuper.ciExpo,by=c("exposure","trans","setting","LT"), all.x = T)



## create table with true HRs --------------------------------------------
temp <- melt(Settings,id.vars = "setting") 
TranslationSettings<-data.table(variable=c('HRSABCur','HRETOPCur','HRSABPrev','HRETOPPrev'),
                                trans=c(2,3,2,3),
                                exposure=c(1,1,2,2))
hr.true<-merge(TranslationSettings,temp,by="variable", all.x = T)
hr.true[,"variable":=NULL]
temp<-expand.grid(trans=1,exposure=c(1,2),setting=1:nsettings,value=1)
hr.true<-rbind(hr.true,temp)
l<-replicate(nLTs, hr.true, simplify = FALSE)
hr.true<-rbindlist(l,idcol="LT")
hr.true[,c("exposure","trans","setting"):=list(factor(exposure),as.character(trans),as.numeric(setting))]

coefs.combined<-merge(coefs.combined,hr.true,by=c("exposure","trans","setting","LT"), all.x = T)



# grouped by setting and LT, both LM and Tdc estimates -------------------------------------------
colsCoefs <- c("landmark"= "#F8766D","tdc"="#1F78B4","true"="black")
pd2<-position_dodge(0.5)
legendp<-matrix(c(0.8,0.9,0.8,0.9,0.8,0.9),nrow=3,byrow = T)
ymin.tr=rep(0.2,3)
ymax.tr=c(2,7,4)

theme_set(theme_minimal(base_size=16))


plotSettings<-function(tr){
  coefs.estsSAB<-coefs.combined[trans==tr]
  for(setti in 1:nsettings){
    coefs.ests<-coefs.estsSAB[setting==setti]
    pl<-ggplot(coefs.ests, aes(x=LM,y=HR.x, color="landmark")) +
      geom_line(aes(x =LM, y = lower.y, color="tdc"), linetype="dashed")+
      geom_line(aes(x =LM, y = upper.y, color="tdc"), linetype="dashed")+
      geom_line(aes(x =LM, y = HR.y, color="tdc"))+
      geom_line(aes(x =LM, y = lowerSuper,  color="landmark" ),linetype="dashed")+
      geom_line(aes(x =LM, y = upperSuper,  color="landmark" ), linetype="dashed")+
      geom_line(aes(x =LM, y = HRSuper,color="landmark"))+
      geom_line(aes(x =LM, y = value,color="true"), linetype="dotted",size=0.9)+
      geom_errorbar(aes(ymax =upper.x, ymin = lower.x, color="landmark"), width=0, position = pd2)+
      geom_point( size=1)+
      theme_minimal()+
      scale_color_manual(name="Model",values=colsCoefs,labels=c("Landmark",
                                                                "TDC","True Value"))+
      xlab("Landmark time")+
      ylab(paste("Hazard ratio for",event.lbls[tr]))+
      ggtitle(paste("Setting",setti))+
      facet_grid(exposure ~ LT, scales="free_y",
                 labeller = labeller(LT=LT.labs,exposure=exposure.labs))+
      theme(strip.background = element_rect(fill="lightgrey"))+
      theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")
    print(pl)
    ggsave(paste0(wdimages,"S",setti+6,"_HR",tr,"Setting",setti,".tiff"), pl,width=6,height=8,dpi=300)
    
  }
}

plotSettings(2) # only for outcome SAB



################################################################################################
## Transition probabilities
probs<-probsCrude
probsCrude.ci<-probs[,list(CIF=mean(CIF),lower=quantile(CIF,probs=0.025),upper=quantile(CIF,probs=0.975)),by=list(LM,exposure,trans,setting,LT)]

probsEtm<-msprobs
probsEtm.ci<-probsEtm[,list(CIF=mean(CIF),lower=quantile(CIF,probs=0.025),upper=quantile(CIF,probs=0.975)),by=list(LM,exposure,trans,setting,LT)]

probs.combined<-rbind(cbind(model=1,probsCrude.ci),cbind(model=2,probsEtm.ci))

## Supermodel #############
probsSuper<-probsSupermodel
probsSuper.ci<-probsSuper[,list(CIF=mean(CIF),lower=quantile(CIF,probs=0.025),upper=quantile(CIF,probs=0.975)),by=list(LM,exposure,trans,setting,LT)]
probsSuper.ci<-cbind(probsSuper.ci,model=1)

probs.combined[,"exposure":=factor(exposure)]
probs.combined<-merge(probs.combined,probsSuper.ci,by=c("model","LM","exposure","trans","setting","LT"),all.x = T)



ymin.tr=c(0.4,0,0)
ymax.tr=c(1,0.4,0.4)
legendp<-matrix(c(0.8,0.2,0.8,0.9,0.8,0.9),nrow=3,byrow = T)

pd2<-position_dodge(0.3)
exposure.labs <- exposure.lbls
names(exposure.labs) <- 0:2

tr<-2 ## only SAB
probs.estsSAB<-probs.combined[trans==tr]
for(setti in 1:nsettings){
    probs.ests<-probs.estsSAB[setting==setti]
    pl<-ggplot(probs.ests, aes(x=LM,y=CIF.x,color=factor(model))) +
      geom_errorbar(aes(ymax =upper.x, ymin = lower.x), width=0, position = pd2)+
      geom_point( size=1.5, position = pd2)+
      theme_minimal()+
      theme(legend.position = "none")+
      xlab("Landmark time")+
      ylab(paste0("Predicted cumulative incidence for ",event.lbls[tr]," at horizon"))+
      ggtitle(paste("Setting",setti))+
      geom_line(aes(x=LM,y=CIF.y))+
      geom_line(aes(x=LM,y=lower.y), linetype="dashed")+
      geom_line(aes(x=LM,y=upper.y), linetype="dashed")+
      geom_line(data=probsTrue[trans==tr&setting==setti],aes(x=LM,y=CIF),color="black", size=0.9, linetype="dotted")+
      facet_grid(LT ~ exposure, scales="free_y",
                 labeller = labeller(LT=LT.labs,exposure=exposure.labs))+
      theme(strip.background = element_rect(fill="lightgrey"))+
      theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")+
      scale_color_hue(name="Model",labels=c("Landmark","Multistate","True Value"))
    
    print(pl)
   ggsave(paste0(wdimages,"S",setti+15,"_Probs","Setting",setti,".tiff"), pl,width=6,height=8,dpi=300)
}  




###############################################
## Tables S: RMSE of the landmark CIFS ########
###############################################

vars<-c("setting","trans","exposure","LT","LM","iteration")

probsCrude.ci<-probsCrude
probsCrude.ci[,"CIFCrude":=CIF]
probsEtm.ci<-msprobs
probsEtm.ci[,"CIFMstate":=CIF]
probsSuper.ci<-probsSupermodel
probsSuper.ci[,"CIFSuper":=CIF]
probsSuper.ci[,"exposure":=as.numeric(as.character(exposure))]

## number of deviations >1
temp<-merge(probsCrude.ci[,c(vars,"CIFCrude"),with=F],probsEtm.ci[,c(vars,"CIFMstate"),with=F],by=vars, all.y = T)
temp2<-merge(temp,probsSuper.ci[,c(vars,"CIFSuper"),with=F],by=vars)
probsTrue[,"setting":=as.numeric(as.character(setting))]
temp3<-setDT(dplyr::left_join(temp2,probsTrue[,c("setting","trans","exposure","LM","CIF"),with=F],
             by=c("setting","trans","exposure","LM"),relationship = "many-to-many"))
temp3[,c("Crude","Super","Mstate"):=list(sqrt(((CIFCrude-CIF)^2))*100,
                             sqrt(((CIFSuper-CIF)^2))*100,
                             sqrt(((CIFMstate-CIF)^2))*100)]
temp3[iteration==670&setting==4&LM==1&LT==1&exposure==2]
temp3[iteration==161&setting==9&LM==1&LT==1&exposure==1]

## absolute error
errorProbsCrude<-abs(temp3$CIFCrude*100-temp3$CIF*100)
errorProbsMstate<-abs(temp3$CIFMstate*100-temp3$CIF*100)
errorProbsSuper<-abs(temp3$CIFSuper*100-temp3$CIF*100)
summary(errorProbsCrude)
summary(errorProbsMstate)
summary(errorProbsSuper)
sum(errorProbsCrude>1,na.rm = T) 
sum(errorProbsMstate>1,na.rm = T)
sum(errorProbsSuper>1,na.rm = T) 

temp3[errorProbsCrude>1]
temp3[errorProbsSuper>1]

RMSE<-temp3[,list(Crude=sqrt(mean((CIFCrude-CIF)^2, na.rm=T))*100,
            Super=sqrt(mean((CIFSuper-CIF)^2, na.rm=T))*100,
            Mstate=sqrt(mean((CIFMstate-CIF)^2, na.rm=T))*100),
      by=list(setting,LT)]
print(xtable(RMSE,digits=c(0,0,0,2,2,2),caption = paste("Root-mean-square error of the estimated landmark CIFs in percent.",
                                                       "LT 1-4 Left-truncation scenarios 1 (Skew normal), 2 (Exp(0.1)), 3 (Exp(0.5)), 4 (Exp(1)).",
                                                       "Crude LM, Crude landmark model. Mstate, multi-state model",
                                                       "Super LM, landmark supermodel."),
                            label = paste("tab: RMSE probs")), include.rownames = F,
      size="\\fontsize{10pt}{10pt}\\selectfont",
      tabular.environment = "longtable", floating = F)


RMSEoutcome<-temp3[,list(Crude=sqrt(mean((CIFCrude-CIF)^2))*100,
                  Super=sqrt(mean((CIFSuper-CIF)^2))*100,
                  Mstate=sqrt(mean((CIFMstate-CIF)^2))*100),
            by=list(trans,setting,LT)]
setorderv(RMSEoutcome,c("trans","LT"))
RMSEoutcome[,c("trans","setting"):=list(event.lbls[trans],factor(setting))]
setnames(RMSEoutcome,old="trans",new="Outcome")
print(xtable(RMSEoutcome,digits=c(0,0,0,0,2,2,2),caption = paste("Root-mean-square error of the estimated landmark CIFs in percent per outcome.",
                                                      "LT 1-4 Left-truncation scenarios 1 (Skew normal), 2 (Exp(0.1)), 3 (Exp(0.5)), 4 (Exp(1)).",
                                                      "Crude LM, Crude landmark model. Mstate, multi-state model",
                                                      "Super LM, landmark supermodel."),
                            label = paste("tab: RMSE probs per outcome")), include.rownames = F,
      size="\\fontsize{10pt}{10pt}\\selectfont",
      tabular.environment = "longtable", floating = F)
