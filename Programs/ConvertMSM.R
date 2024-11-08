

###################################################################################
## Code to transform the imputet data into the required multistate model format ###
###################################################################################

# inputs: data.table in wide format (dwide) with columns
# "idChild": identifying the fetus
# "entry","exit": study entry and exit
# "from_1", "to_1", "from_2", "to_2": times of first and second exposure (NA if no exposure occured),
# covariates (covs)
# integer indicating whether or not to combine the currently and previously exposed state (model=2,1)
# ouput: dataset in multisate format

fconvert_to_MSM<-function(dwide,covs,model=1){
  d<-setDT(dwide)
  
  # first and second exposure are directly back to back-> combine first and second exposre
  d[to_1>=from_2]
  d[to_1>=from_2,`:=`("to_1"=to_2, "from_2"=NA,"to_2"=NA)] 
  
  dlong_list<-lapply(1:nrow(d),ftransitions_wide_to_long,d=d,covs=covs,model=model)

  dlong<-rbindlist(dlong_list)
  setnames(dlong, old=c("entry","exit","Tstart","Tstop"), new=c("study.entry", "study.exit","entry","exit"))
  
  # only observed transitions - excluding the ones that occurred before study entry
  dMSM<-subset(dlong,exit>study.entry | exit==study.exit)
  dMSM<-cbind("number"=c(1:nrow(dMSM)),dMSM)
  
  #temp=subset(dMSM, entry<study.entry & exit!=study.exit) # states with study study.entry between to and from
  temp<-subset(dMSM, entry<study.entry) # transitions with exposure start before study entry (but exposure end after study entry)
  temp$entry<-temp$study.entry
  
  temp2<-dMSM[number%in%setdiff(dMSM$number,temp$number)]
  dMSM<-rbind(temp,temp2)
  dMSM<-dMSM[order(dMSM$number),]
  dMSM$number<-NULL
  
  dMSM<-dMSM[,c("idChild","idWoman","study.entry","study.exit", "entry","exit","from","to",covs),with=F]
  
  return(dMSM)
}



ftransitions_wide_to_long<-function(i,d,covs,model=1){
  x<-d[i,]
  
  curstate<-0 # never exposed before at study entry
  tcur<-0
  
  
  df<-NULL # data frame with transition states and times for women i
  if(!is.na(x$from_1)){ 
    nextstate<-1 # currently exposed
    nextt<-x$from_1
    df<-rbind(df,c("Tstart"=tcur,"Tstop"=nextt,"from"=curstate,"to"=nextstate))
    tcur<-x$from_1
    curstate<-nextstate
  }
  
  if(model==1){  # differentiate currently and previously exposed state
    if(!is.na(x$to_1)  ){
      if(x$to_1>=x$exit){ # exposure ends after study exit
        nextstate<-x$outcome+2
        nextt<-x$exit # cut exposure duration at study exit
      }else{
        nextstate<-2
        nextt<-x$to_1
      }
      df<-rbind(df,c("Tstart"=tcur,"Tstop"=nextt,"from"=curstate,"to"=nextstate))
      tcur<-nextt
      curstate<-nextstate
    }
    
    if(!is.na(x$from_2)){ # second exposure
      nextstate<-1
      nextt<-x$from_2
      df<-rbind(df,c("Tstart"=tcur,"Tstop"=nextt,"from"=curstate,"to"=nextstate))
      tcur<-nextt
      curstate<-nextstate
    }
    
    if(!is.na(x$to_2)){
      if(x$to_2>=x$exit){
        nextstate<-x$outcome+2
        nextt<-x$exit
      }else{
        nextstate<-2
        nextt<-x$to_2
      }
      df<-rbind(df,c("Tstart"=tcur,"Tstop"=nextt,"from"=curstate,"to"=nextstate))
      tcur<-nextt
      curstate<-nextstate
    }
  }
  
  
  
  if(tcur<x$exit){ # no exposure
    xx<-data.frame("Tstart"=tcur,"Tstop"=x$exit, "from"=curstate,"to"=x$outcome+2)
    df<-rbind(df,xx)
  }
  df
  
  df2<-cbind(df,x[,c("entry","exit","idChild","idWoman",covs),with=F])
  
  return(df2)
}

