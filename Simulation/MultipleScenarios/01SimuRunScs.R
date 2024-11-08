#####################################################################################################
# Simulation of data and estimation of hazard and cumulative incidences in additional settings ######
#####################################################################################################

# changes on the cluster: pacioli<-T, dopar instead of do, nit, n

library(data.table)
library(mstate)
library(doParallel)

### needs 2 directories under wd: RunFiles (load cumHazSimu.Rda, LMSimulation.R, AnaLandmark.R), Results


n<-6000 # number of women per iteration
nparallel<-1

pacioli<-T # indicating whether the code is run on a cluster (pacioli) or a local test version is run for spot checks
if(pacioli==T){
  wd<-getwd()
  registerDoParallel(cores=4)
  nit<-1000 # number of iterations per run
}
# if(pacioli==F){ # local test version for spot checks
#   wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
#   # reduced number of iterations for running the simulation on a local computer
#   nit<-2
# }

runfiles.path<-paste0(wd,"/RunFiles/") # to load cumHazSimu.Rda and for LMSimulation.R (which calls AnaLandmark.R)
Settings<-readRDS(paste0(runfiles.path,"Settings.Rds"))
source(paste0(runfiles.path,"LMSimulationScs.R"),local=T)

nLTdists<-length(LTpars) # number of left truncation distributions
if(pacioli==F) Settings<-Settings[1:2,] # reduced number of settings for spot checks
settings<-unique(Settings$setting)
RNGkind(sample.kind = "Rounding")  # match random number generation prior to R 3.6.0



reslist<-foreach(iteration= seq_len(nit),.combine=list) %do%{
  reslistSet<-foreach(sc=iter(settings),.combine=list) %do%{
    ParsSetting<-Settings[setting==sc,]
    ##########################################################################################################
    ## Simulation ############################################################################################
    ## save results
    out1<-simMSM(iteration,ParsSetting,nit,nparallel)
    
    #---add LT times --------------------------------------------
    reslistLT<-foreach(LTi= seq_len(nLTdists)) %do%{
      simData<-appendScenarios(out1,LTi, iteration,nit,nparallel)
      msData<-copy(simData)
      simData<-simData[status==1]
      simDataTranslist<-with(simData,table(from,to))
      simDataTranslist<-as.data.frame.matrix(simDataTranslist)
      simDataTranslist$from<-rownames(simDataTranslist)
      ##########################################################################################################
      ## TDC analysis ##########################################################################################
      EstsTdclist<-lapply(4:6, estsTdc,TdcData=simData)
      EstsTdc<-rbindlist(EstsTdclist,idcol="trans")
      
      ##########################################################################################################
      ## Landmark analysis #####################################################################################
      LMlong<-prepare.LMData(simData)
      
      # save for description
      lmDataTrans<-lmCounts(LMlong)
      
      
      ## Estimation #############
      ## Separate ----------------------------------------------------------------------------------
      EstsSeparate.list<-lmEstsSeparate(LMlong)
      coefsCrude<-EstsSeparate.list[[1]]
      probsCrude<-EstsSeparate.list[[2]]
      
      
      ## Supermodel ----------------------------------------------------------------------------------
      EstsSupermodel.list<-lmEstsSupermodel(LMlong)
      
      coefsSupermodelAll<-EstsSupermodel.list[[1]]
      coefsSupermodel<-EstsSupermodel.list[[2]]
      probsSupermodel<-EstsSupermodel.list[[3]]
      
      
      #############################################################################################
      ### Non-parametric estimates ##################################################################
      #############################################################################################
      msprobs<-msEstimates(msData)
      
      
      ##############################################################################################
      reslistLT<-list(simDataTranslist=simDataTranslist,
                      lmDataTrans=lmDataTrans,
                      coefsTdc=EstsTdc,
                      coefsCrude=coefsCrude,
                      probsCrude=probsCrude,
                      coefsSupermodelAll=coefsSupermodelAll, 
                      coefsSupermodel=coefsSupermodel,
                      probsSupermodel=probsSupermodel,
                      msprobs=msprobs)
    }
    
    simDataTranslisti<-lapply(reslistLT,function(l)l$simDataTranslist)
    lmDataTranslisti<-lapply(reslistLT,function(l)l$lmDataTrans)
    msprobslisti<-lapply(reslistLT,function(l)l$msprobs)
    coefsTdclisti<-lapply(reslistLT,function(l)l$coefsTdc)
    coefsCrudelisti<-lapply(reslistLT,function(l)l$coefsCrude)
    probsCrudelisti<-lapply(reslistLT,function(l)l$probsCrude)
    coefsSupermodelAlllisti<-lapply(reslistLT,function(l)l$coefsSupermodelAll)
    coefsSupermodellisti<-lapply(reslistLT,function(l)l$coefsSupermodel)
    probsSupermodellisti<-lapply(reslistLT,function(l)l$probsSupermodel)
    
    simDataTranslisti<-rbindlist(simDataTranslisti,idcol = "LT")
    lmDataTranslisti<-rbindlist(lmDataTranslisti,idcol = "LT")
    msprobslisti<-rbindlist(msprobslisti, idcol = "LT")
    coefsTdclisti<-rbindlist(coefsTdclisti,idcol = "LT")
    coefsCrudelisti<-rbindlist(coefsCrudelisti,idcol = "LT")
    probsCrudelisti<-rbindlist(probsCrudelisti,idcol = "LT")
    coefsSupermodelAlllisti<-rbindlist(coefsSupermodelAlllisti,idcol = "LT")
    coefsSupermodellisti<-rbindlist(coefsSupermodellisti,idcol = "LT")
    probsSupermodellisti<-rbindlist(probsSupermodellisti,idcol = "LT")
    
    reslist_LT<-list(simDataTranslist=simDataTranslisti,
                     lmDataTrans=lmDataTranslisti,
                     msprobs=msprobslisti,
                     coefsTdc=coefsTdclisti,
                     coefsCrude=coefsCrudelisti,
                     probsCrude=probsCrudelisti,
                     coefsSupermodelAll=coefsSupermodelAlllisti,
                     coefsSupermodel=coefsSupermodellisti,
                     probsSupermodel=probsSupermodellisti)  
    return(reslist_LT)
    
  }
  
  ##########################################################################################################
  simDataTranslistSet<-lapply(reslistSet,function(l)l$simDataTranslist)
  lmDataTranslistSet<-lapply(reslistSet,function(l)l$lmDataTrans)
  msprobslistSet<-lapply(reslistSet,function(l)l$msprobs)
  coefsTdclistSet<-lapply(reslistSet,function(l)l$coefsTdc)
  coefsCrudelistSet<-lapply(reslistSet,function(l)l$coefsCrude)
  probsCrudelistSet<-lapply(reslistSet,function(l)l$probsCrude)
  coefsSupermodelAlllistSet<-lapply(reslistSet,function(l)l$coefsSupermodelAll)
  coefsSupermodellistSet<-lapply(reslistSet,function(l)l$coefsSupermodel)
  probsSupermodellistSet<-lapply(reslistSet,function(l)l$probsSupermodel)
  
  simDataTranslistSet<-cbind("iteration"=iteration,rbindlist(simDataTranslistSet,idcol = "setting"))
  lmDataTranslistSet<-cbind("iteration"=iteration,rbindlist(lmDataTranslistSet,idcol = "setting"))
  msprobslistSet<-cbind("iteration"=iteration,rbindlist(msprobslistSet, idcol = "setting"))
  coefsTdclistSet<-cbind("iteration"=iteration,rbindlist(coefsTdclistSet,idcol = "setting"))
  coefsCrudelistSet<-cbind("iteration"=iteration,rbindlist(coefsCrudelistSet,idcol = "setting"))
  probsCrudelistSet<-cbind("iteration"=iteration,rbindlist(probsCrudelistSet,idcol = "setting"))
  coefsSupermodelAlllistSet<-cbind("iteration"=iteration,rbindlist(coefsSupermodelAlllistSet,idcol = "setting"))
  coefsSupermodellistSet<-cbind("iteration"=iteration,rbindlist(coefsSupermodellistSet,idcol = "setting"))
  probsSupermodellistSet<-cbind("iteration"=iteration,rbindlist(probsSupermodellistSet,idcol = "setting"))
  
  reslist_iteration<-list(simDataTrans=simDataTranslistSet,
                lmDataTrans=lmDataTranslistSet,
                msprobs=msprobslistSet,
                coefsTdc=coefsTdclistSet,
                coefsCrude=coefsCrudelistSet,
                probsCrude=probsCrudelistSet,
                coefsSupermodelAll=coefsSupermodelAlllistSet,
                coefsSupermodel=coefsSupermodellistSet,
                probsSupermodel=probsSupermodellistSet)
  return(reslist_iteration)
}


simDataTrans<-rbindlist(lapply(reslist,function(l)l$simDataTrans))
lmDataTrans<-rbindlist(lapply(reslist,function(l)l$lmDataTrans))
msprobs<-rbindlist(lapply(reslist,function(l)l$msprobs))
coefsTdc<-rbindlist(lapply(reslist,function(l)l$coefsTdc))
coefsCrude<-rbindlist(lapply(reslist,function(l)l$coefsCrude))
probsCrude<-rbindlist(lapply(reslist,function(l)l$probsCrude))
coefsSupermodelAll<-rbindlist(lapply(reslist,function(l)l$coefsSupermodelAll))
coefsSupermodel<-rbindlist(lapply(reslist,function(l)l$coefsSupermodel))
probsSupermodel<-rbindlist(lapply(reslist,function(l)l$probsSupermodel))

# save estimation results
results.path<-ifelse(pacioli==T,paste0(wd,"/Results/Cluster/"),paste0(wd,"/Results/Test/"))
saveRDS(simDataTrans,file=paste0(results.path,"simDataTrans.Rds"))
saveRDS(lmDataTrans,file=paste0(results.path,"lmDataTrans.Rds"))
saveRDS(msprobs,file=paste0(results.path,"msprobs.Rds"))
saveRDS(coefsTdc,file=paste0(results.path,"coefsTdc.Rds"))
saveRDS(coefsCrude,file=paste0(results.path,"coefsCrude.Rds"))
saveRDS(probsCrude,file=paste0(results.path,"probsCrude.Rds"))
saveRDS(coefsSupermodelAll,file=paste0(results.path,"coefsSupermodelAll.Rds"))
saveRDS(coefsSupermodel,file=paste0(results.path,"coefsSupermodel.Rds"))
saveRDS(probsSupermodel,file=paste0(results.path,"probsSupermodel.Rds"))







