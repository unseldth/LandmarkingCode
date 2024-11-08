## Load packages and functions
library(data.table)
library(mstate)
library(xtable)
library(ggplot2)
library(latex2exp)
library(scales)
library(RColorBrewer)
library(doParallel)

data.path<-paste0(wd,"/Data/") 
wdimages<-paste0(wd,"/Figures/") 
wdresults<-paste0(wd,"/Results/") 
programs.path<-paste0(wd,"/Programs/")

event.lbls<-c("Live birth","SAB", "ETOP", "Stillbirth")

expo.lbls<-c("No prior exposure","Prior exposure")
names(expo.lbls)<-0:1

age.lbls<-c("20-35 years","<=20 years", ">35 years" ) 
age.labs<-paste("Age",c("moderate","low","high"))
names(age.labs)<-0:2


trans.labs<-event.lbls
names(trans.labs)<-1:3
models.lbls<-c("Unadjusted","Adjusted")
names(models.lbls)<-1:2
probs.lbls<-c("Unadjusted","Moderate age group")
names(probs.lbls)<-1:2

base.covs<-c("age.gr","prev_preg",
             "prev_para","prev_fet_loss","prev_child_bd")

state.names.long<-c(expo.lbls,event.lbls)
nstates<-length(state.names.long)
state.names<-as.character(c(0,1,3:nstates))
names(state.names)

tra.num<-transMat(x = list(c(2, 3:nstates), c(3:nstates),c(), c(), c(),c()), names = state.names)


tra<-matrix(F,nstates,nstates)
tra[1,c(2:nstates)]<-T
tra[2,c(3:nstates)]<-T
colnames(tra)<-rownames(tra)<-state.names


nevents<-length(event.lbls)
tmatCompete <- trans.comprisk(nevents, names = event.lbls)
