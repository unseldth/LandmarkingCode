# Set path for input and output to current working directory
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)

# Load global variables and helper functions
source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"HelpFunctionsAnaGyrase.R"), local = T)
source(paste0(programs.path,"ConvertMSM.R"), local = T)

imp_anonymized<-readRDS(file=paste0(data.path,"imp_anonymized.Rds"))


################################################################################
## Convert to multistate model and save results in data folder #################
################################################################################
head(imp_anonymized)


# distinct currently and previously exposed state
dMSM<-fconvert_to_MSM(imp_anonymized,covs=c("age.gr"))
setorder(dMSM,idWoman,idChild)

# combined currently and previously exposed state
dMSM2<-fconvert_to_MSM(imp_anonymized,covs=c("age.gr"),model = 2) 

table(dMSM$from,dMSM$to)
table(dMSM2$from,dMSM2$to)


saveRDS(dMSM,file=paste0(data.path,"dMSM_anonymized.Rds"))
saveRDS(dMSM2,file=paste0(data.path,"dMSM2.Rds"))


##########################################################################################################
## Transform dMSM2 data in the format needed for mstate with status=0,1 and trans column #################
##########################################################################################################

# function adapted from the msprep() function in the mstate package
# input: non-absorbing state, output: binary vector indicating the tansition status
mymsprep<-function(transition){
  s<-subset(dMSM2,from==transition)
  possible.trans<-names(which(tra[transition,]==T))
  # append covariates and number by copying info in s
  s2<-s[rep(seq_len(nrow(s)), each = length(possible.trans)), ] 
  s2$to2<-rep(possible.trans,times=nrow(s))
  to2<-rep(possible.trans,times=nrow(s))
  y<-as.character(s2$to)
  s2$status<-ifelse(s2$to2==y,1,0)
  return(s2)
}
non.absorbin<-names(which(rowSums(tra)>0))
l<-lapply(non.absorbin,mymsprep)
dMSMexpand<-do.call("rbind",l)
dMSMexpand$to<-NULL
setnames(dMSMexpand,old="to2",new="to")

dMSMexpand$trans<-tra.num[cbind(as.character(dMSMexpand$from),as.character(dMSMexpand$to))]
dMSMexpand[,`:=`(to=factor(to),trans=factor(trans))]
dMSMexpand<-dMSMexpand[,c( "idChild","idWoman",
                           "entry" ,"exit" ,"from","to","status","trans",
                           "age.gr"),with=F]
dMSMexpand<-dMSMexpand[order(dMSMexpand$idChild),]

saveRDS(dMSMexpand,file=paste0(data.path,"dMSMexpand.Rds"))


###########################################################################################
## transform dMSM2 data in the format needed for landmarking with msfit() #################
###########################################################################################


temp<-copy(imp_anonymized)
temp[,"to":=outcome+2]
# temp<-temp[to!=6] # exclude stillbirths
head(temp)

LMlong<-prepare.LMData(temp, covs=c("age.gr"), idOrig = "idChild")
table(LMlong$LM)
# table(LMlong%>%filter(status==1)%>%select(to))


LMlong<-LMlong[with(LMlong, order(idChild,Tstart,trans)), ]
idcols<-c("idChild","idWoman")
LMlong<-LMlong[,c(idcols,setdiff(colnames(LMlong),idcols))]
head(LMlong)


attr(LMlong, "trans") <- tmatCompete
saveRDS(LMlong,paste0(data.path,"/LMlong.Rds"))

