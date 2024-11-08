library(data.table)
library(mstate)


wdimages<-paste0(wd,"/../Figures/Simulation/")  
runfiles.path<-paste0(wd,"/RunFiles/") # to load cumHazSimu.Rds and for LMSimulation.R (which calls AnaLandmark.R)
data.path<-paste0(wd,"/Data/")  
results.path<-paste0(wd,"/Results/Cluster/")  # saved estimation results



## Variable labels ####
state.lbls<-paste("State",1:6-1)

state.lbls<-state.lbls2<-paste("State",1:6-1)
event.lbls<-state.lbls[4:6]
exposure.lbls<-state.lbls[1:3]
hr.lbls<-paste("Exposure",c(1:2)," vs. 0")
hr.labs<-hr.lbls
names(hr.labs)<-1:2