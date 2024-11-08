##############################################################################
# Simulation of data and estimation of hazard and cumulative incidences ######
##############################################################################

# changes on the cluster: pacioli<-T, dopar instead of do, nparallel, nitPerRun, n
library(data.table)
library(mstate)
library(doParallel)
### needs 2 directories under wd: RunFiles (load cumHazSimu.Rds, LMSimulation.R, AnaLandmark.R), Results


nparallel<-1 # id of parallel run, change to 1,2,...,10 when running parallel jobs on the cluster
n<-6000 # number of women per iteration
scenarios<-1:3
scenariosSeperate<-1:2

pacioli<-T # indicating whether the code is run on a cluster (pacioli) or a reduced local test version for spot checks
if(pacioli==T){
  wd<-getwd()
  registerDoParallel(4)
  Nparallel<-10 # number of parallel runs
  nitPerRun<-1000
  nit<-Nparallel*nitPerRun
}
# if(pacioli==F){ # local test version for spot checks
#   wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
#   # reduced number of iterations for running the simulation on a local computer
#   Nparallel<-1
#   nitPerRun<-2
#   nit<-Nparallel*nitPerRun
# }

runfiles.path<-paste0(wd,"/RunFiles/") # to load cumHazSimu.Rds and for LMSimulation.R (which calls AnaLandmark.R)
source(paste0(runfiles.path,"LMSimulation.R"),local=T)
RNGkind(sample.kind = "Rounding")  # match random number generation prior to R 3.6.0

reslist<-foreach(iteration= seq_len(nitPerRun),.combine=list) %do%{
  ##########################################################################################################
  ## Simulation ############################################################################################
  out1<-simMSM(iteration,nitPerRun,nparallel)
  out1[Tstop==Inf,c("Tstop","status"):=list(45,ifelse(to==4,1,0))]
  out1[,"duration":=Tstop-Tstart]
  
  # add scenarios 2 and 3
  simData<-appendScenarios(out1)
  msData<-copy(simData)
  simData<-simData[status==1]
  
  ## save results
  simDataTranslist<-lapply(1:length(scenarios),function(sc) with(subset(simData,scenario==sc),table(from,to)))
  
  
  ##########################################################################################################
  ## Landmark analysis #####################################################################################
  LMlong<-prepare.LMData(data.MSM=simData)
  
  #save for description
  lmDataTranslist<-lapply(1:length(scenarios), lmCounts,LMlong=LMlong)
  lmDataTrans<-rbindlist(lmDataTranslist,idcol="scenario")


  ## Estimation #############
  ## Crude ----------------------------------------------------------------------------------
  EstsSeparate.list<-lmEstsSeparate(LMlong)

  coefsCrude<-EstsSeparate.list[[1]]
  probsCrude<-EstsSeparate.list[[2]]


  
  ## Supermodel ----------------------------------------------------------------------------------
  EstsSupermodel.list<-lmEstsSupermodel(LMlong)
  
  coefsSupermodelAll<-EstsSupermodel.list[[1]]
  coefsSupermodel<-EstsSupermodel.list[[2]]
  probsSupermodel<-EstsSupermodel.list[[3]]
  
  
  ############################################################################################
  ## Non-parametric estimates ##################################################################
  ############################################################################################
  msprobslist<-lapply(scenarios,msEstimates,msData=msData)
  msprobs<-rbindlist(msprobslist, idcol = "scenario")

  
  res<-list(simDataTranslist=simDataTranslist,
       lmDataTrans=lmDataTrans,
       probsCrude=probsCrude,
       coefsCrude=coefsCrude,
       coefsSupermodelAll=coefsSupermodelAll,
       coefsSupermodel=coefsSupermodel,
       probsSupermodel=probsSupermodel,
       msprobs=msprobs)

  return(res)
}


if(pacioli) saveRDS(reslist,paste0(wd,"/Results/Cluster/reslist",nparallel,".Rds"))
if(pacioli==F) saveRDS(reslist,paste0(wd,"/Results/Test/reslist",nparallel,".Rds"))





