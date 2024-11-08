###################################################################
# Representation of the hazards by parametric functions ##########
###################################################################


###Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
source(paste0(wd,"/GlobalVars.R"))
runfiles.path<-paste0(wd,"/Simulation/RunFiles/")


library(mvna)
library(sn)
library(latex2exp)
library(flexsurv)
library(scales)
library(RColorBrewer)
library(dynpred)
library(gridExtra)



dMSM<-readRDS(file=paste0(wd,"/../Data/dMSM_anonymized.Rds"))




clrs<-c(hue_pal()(3),"indianred1","plum") 
clrs.paired<-brewer.pal(3, "Paired")

# state.lbls<-c("Never exposed before","Currently exposed","Previously exposed","Life Birth","SAB","ETOP")

state.names<-as.character(0:5)
state.coding<-cbind(state.names,state.lbls)

###Time grid
t<-seq(0.00,50,0.05)
cutt<-45 # no events after this week

#################################################################################
## Model for left-truncation times ###############################################
#################################################################################
ordered_data <- as.data.frame(dMSM[order(dMSM$idChild, dMSM$entry),])
data.first.state <- ordered_data[!duplicated(ordered_data$idChild),]

LTpars=list(list(distr="SN",xi=9,omega = 4.3,alpha = -1),
            list(distr="Exp",rate=0.1),
            list(distr="Exp",rate=0.5),
            list(distr="Exp",rate=1.5))
parsL<-LTpars[[1]]
# density of observed study entries
hist(data.first.state$study.entry,freq = F, main=NULL,xlab="Gestational week")

curve(dsn(x,xi=parsL$xi,omega = parsL$omega,alpha = parsL$omega),add=TRUE,lwd=2,col="orange")
for(i in 2:4){
   LTdist<-LTpars[[i]]
   curve(dexp(x,rate=LTdist$rate),add=TRUE,lwd=2,col=i+1)
}


# empirical distribution function of observed study entries
FL<-ecdf(data.first.state$study.entry)
plot(FL,xlim=c(0,cutt),xlab="Gestational week",ylab="Distribution function (LT times)",lwd=2,main="")
curve(psn(x,xi=parsL$xi,omega = parsL$omega,alpha = parsL$omega),xlim = c(-10,cutt),add=TRUE,lwd=2,col="orange",xlab="Gestational week")
for(i in 2:4){
   LTdist<-LTpars[[i]]
   curve(pexp(x,rate=LTdist$rate),add=TRUE,lwd=2,col=i+1)
}



#################################################################################
## Models for transition hazards ################################################
#################################################################################

## NAEs from gyrase data
tmat<-!is.na(transMat(x = list(c(1, 3:5)+1, c(2,3:5)+1,c(1,3:5)+1, c(), c(), c()), 
                      names = 0:5))




#####################################################################################
## Transitions to the outcomes  ####################################################
## Life birth ----------------------------------------------------------------------
pars.life<-data.frame(shape= c(0.45,0.45,0.45),rate=rep(5*10^(-9),3))

A03<-sapply(pmin(t,cutt),Hgompertz, shape=pars.life[1,1],rate=pars.life[1,2])
A13<-sapply(pmin(t,cutt),Hgompertz, shape=pars.life[2,1],rate=pars.life[2,2])
A23<-sapply(pmin(t,cutt),Hgompertz, shape=pars.life[3,1],rate=pars.life[3,2])



## hazards
hgompertz2<-function(t,shape,rate)hgompertz(t,shape,rate)*1*(t<=cutt)
ymax=3
a03<-sapply(pmin(t,cutt), hgompertz2,shape=pars.life[1,1],rate=pars.life[1,2])
a13<-sapply(pmin(t,cutt), hgompertz2,shape=pars.life[2,1],rate=pars.life[2,2])
a23<-sapply(pmin(t,cutt), hgompertz2,shape=pars.life[3,1],rate=pars.life[3,2])

## SAB ----------------------------------------------------------------------
weibull.cumhaz<-function(t,shape,scale){scale*t^shape}
pars.SAB<-data.frame(shape=c(0.45,0.45,0.45),scale=c(0.03,0.15,0.10))

cuttSAB<-23
A04<-sapply(pmin(t,cuttSAB),weibull.cumhaz, shape=pars.SAB[1,1],scale=pars.SAB[1,2])
A14<-sapply(pmin(t,cuttSAB),weibull.cumhaz, shape=pars.SAB[2,1],scale=pars.SAB[2,2])
A24<-sapply(pmin(t,cuttSAB),weibull.cumhaz, shape=pars.SAB[3,1],scale=pars.SAB[3,2])

## hazards
alpha.weibull<-function(t,shape,scale){shape*scale*t^(shape-1)*1*(t<=cuttSAB)}

a04<-sapply(pmin(c(0.000001,t[-1]),cuttSAB), alpha.weibull,shape=pars.SAB[1,1],scale=pars.SAB[1,2])
a14<-sapply(pmin(c(0.000001,t[-1]),cuttSAB), alpha.weibull,shape=pars.SAB[2,1],scale=pars.SAB[2,2])
a24<-sapply(pmin(c(0.000001,t[-1]),cuttSAB), alpha.weibull,shape=pars.SAB[3,1],scale=pars.SAB[3,2])


## ETOP ----------------------------------------------------------------------
HexpETOP<-function(t,steps,t.limits) {
   #steps=ta;t.limits=tm
   steps[1]*(min(t,t.limits[2])-t.limits[1])*1*(t>=t.limits[1])+
      steps[2]*(min(t,t.limits[3])-t.limits[2])*1*(t>t.limits[2])+
      steps[3]*(min(t,t.limits[4])-t.limits[3])*1*(t>t.limits[3])+
      steps[4]*(min(t,t.limits[5])-t.limits[4])*1*(t>t.limits[4])+
      steps[5]*(min(t,t.limits[6])-t.limits[5])*1*(t>t.limits[5])+
      steps[6]*(min(t,t.limits[7])-t.limits[6])*1*(t>t.limits[6])
}
tm <- c(0,2,5,12,23,30,cutt)
ta <- c(1/1000,1/100,4/100,1/1000,1/10000,0)/2
pars.ETOP<-matrix(c(ta,2*ta,1.5*ta),nrow=3,byrow=T)

A05<-sapply(pmin(t,cutt),  HexpETOP,steps=pars.ETOP[1,],t.limits=tm)
A15<-sapply(pmin(t,cutt),  HexpETOP,steps=pars.ETOP[2,],t.limits=tm)
A25<-sapply(pmin(t,cutt),  HexpETOP,steps=pars.ETOP[3,],t.limits=tm)

## hazards
a05<-evalstep(time=tm[-length(tm)], stepf=pars.ETOP[1,], newtime=t, subst=0)
a15<-evalstep(time=tm[-length(tm)], stepf=pars.ETOP[2,], newtime=t, subst=0)
a25<-evalstep(time=tm[-length(tm)], stepf=pars.ETOP[3,], newtime=t, subst=0)



###################################################################################################
## Transitions between the exposure states #######################################################

rate.never<-0.02
rate.prev<-0.01 # a bit lower rate because having to take the fluoroquinolone again after just having taken it is rarer
rate.expo<-1

A01<-sapply(pmin(t,cutt),Hexp, rate=rate.never)
A21<-sapply(pmin(t,cutt),Hexp, rate=rate.prev)
A12<-sapply(pmin(t,cutt),Hexp, rate=rate.expo)


###################################################################################################
## Create matrix of cumulative hazards for mssample ################################################

cum.haz<-rbind(data.frame(time=t,Haz=A01,trans=as.integer(1)),
               data.frame(time=t,Haz=A03,trans=as.integer(2)),
               data.frame(time=t,Haz=A04,trans=as.integer(3)),
               data.frame(time=t,Haz=A05,trans=as.integer(4)),
               data.frame(time=t,Haz=A12,trans=as.integer(5)),
               data.frame(time=t,Haz=A13,trans=as.integer(6)),
               data.frame(time=t,Haz=A14,trans=as.integer(7)),
               data.frame(time=t,Haz=A15,trans=as.integer(8)),
               data.frame(time=t,Haz=A21,trans=as.integer(9)),
               data.frame(time=t,Haz=A23,trans=as.integer(10)),
               data.frame(time=t,Haz=A24,trans=as.integer(11)),
               data.frame(time=t,Haz=A25,trans=as.integer(12)))
saveRDS(cum.haz,file=paste0(wd,"/RunFiles/cumHazSimu.Rds"))





###################################################################################################
### Figure S2 Comparsion of hazards to all possible destination states from each exposure state  ##
###################################################################################################
cum.haz<-data.table(rbind(data.frame(time=t,Haz=A03,from=0,to=3,type="Haz" ),
                          data.frame(time=t,Haz=A13,from=1,to=3,type="Haz" ),
                          data.frame(time=t,Haz=A23,from=2,to=3,type="Haz" ),
                          data.frame(time=t,Haz=A04,from=0,to=4,type="Haz" ),
                          data.frame(time=t,Haz=A14,from=1,to=4,type="Haz" ),
                          data.frame(time=t,Haz=A24,from=2,to=4,type="Haz" ),
                          data.frame(time=t,Haz=A05,from=0,to=5,type="Haz" ),
                          data.frame(time=t,Haz=A15,from=1,to=5,type="Haz" ),
                          data.frame(time=t,Haz=A25,from=2,to=5,type="Haz" ),
                          
                          data.frame(time=t,Haz=A01,from=0,to=1,type="Haz" ),
                          data.frame(time=t,Haz=A12,from=1,to=2,type="Haz" ),
                          data.frame(time=t,Haz=A21,from=2,to=1,type="Haz" )))


# from.labs <- paste("From",c("never exposed before","currently exposed","previously exposed"))
from.labs<-paste("Transition from State",0:2)
names(from.labs) <- c(0:2)


pl<-ggplot(cum.haz,aes(x=time,y = Haz,color=factor(to)))+
   geom_line()+
   scale_y_sqrt()+
   scale_color_hue(name="Transition to", labels=state.lbls[-1])+
   ylab("Cumulative Hazard")+
   xlab("Gestational week")+
   theme(strip.background = element_rect(fill="lightgrey"))+
   facet_grid(from~., scales = "free_y",labeller = labeller(from=from.labs))
pl
ggsave(paste0(wdimages,"S2_TransfromComparison.tiff"), pl,width=8,height=8,dpi=300)



#################################################################################
## Models for transition hazards in additional settings #########################
#################################################################################

cum.haz<-rbind(data.frame(time=t,Haz=A01,trans=as.integer(1)),
               data.frame(time=t,Haz=A03,trans=as.integer(2)),
               #data.frame(time=t,Haz=A04,trans=as.integer(3)),
               #data.frame(time=t,Haz=A05,trans=as.integer(4)),
               #data.frame(time=t,Haz=A12,trans=as.integer(5)),
               data.frame(time=t,Haz=A13,trans=as.integer(6)),
               #data.frame(time=t,Haz=A14,trans=as.integer(7)),
               #data.frame(time=t,Haz=A15,trans=as.integer(8)),
               data.frame(time=t,Haz=A21,trans=as.integer(9)),
               data.frame(time=t,Haz=A23,trans=as.integer(10))
               #data.frame(time=t,Haz=A24,trans=as.integer(11)),
               #data.frame(time=t,Haz=A25,trans=as.integer(12))
)
saveRDS(cum.haz,file=paste0(wd,"/MultipleScenarios/RunFiles/cumHazSimuScs.Rds"))



#######################################################################################
### Parameters for the additional settings -------------------------------------------
HRSABCur=c(2,5)
HRETOPCur=2
RateCurPrev=c(0.1,0.5,1)
Settings<-data.table(expand.grid(HRSABCur=HRSABCur,
                                 HRETOPCur=HRETOPCur,
                                 RateCurPrev=RateCurPrev))
Settings<-rbind(Settings,data.table(HRSABCur=5,RateCurPrev=1,HRETOPCur=1)) # neutralize HR ETOP
Settings<-Settings[order(HRSABCur,RateCurPrev)]

Settings[,c("HRSABPrev","HRETOPPrev"):=list(
   HRSABCur*2/3,
   HRETOPCur)]

# same HR for current and previous exposure
Settings<-rbind(Settings,data.table(HRSABCur=5,RateCurPrev=1,HRETOPCur=2,HRSABPrev=5,HRETOPPrev=2))

Settings[,"setting":=1:nrow(Settings)]
Settings<-Settings[,c("setting",'HRSABCur','HRETOPCur','RateCurPrev','HRSABPrev','HRETOPPrev')]
Settings[,"pStartNever":=0.6]
Settings<-rbind(Settings,Settings[setting==6])
Settings[9,c("setting","pStartNever"):=list(9,0.8)]
Settings[,"setting":=factor(setting)]

saveRDS(Settings,file=paste0(wd,"/MultipleScenarios/RunFiles/Settings.Rds"))


