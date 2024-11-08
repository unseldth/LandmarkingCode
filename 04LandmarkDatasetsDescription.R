## Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"HelpFunctionsAnaGyrase.R"), local = T)

## load packages and functions
library(compareGroups)
library(dplyr)
library(gridExtra)

theme_set(theme_bw(base_size=16))

dMSM2<-readRDS(paste0(data.path,"dMSM2.Rds"))
LMlong<-readRDS(paste0(data.path,"/LMlong.Rds"))

clrs2<-brewer.pal(n = 6,"BrBG")[c(2,1,3:6)]
clrs.paired<-brewer.pal(3, "Paired")




##############################
## Set-up landmark dataset ###
##############################

data<-data.table(LMlong)
data<-data[status==1]


############################################################
## Figure 2: General description of landmark datasets ######
############################################################

nrisk<-data[,list("Never"=sum(exposure=="0"),"Before"=sum(exposure=="1")), by=LM]
nrisk.m <- melt(nrisk,id.vars = "LM") 
# nrisk.m=nrisk.m[order(variable, decreasing = T)]

nevent<-data[,list("Live Birth"=sum(trans==1),"SAB"=sum(trans==2), "ETOP"=sum(trans==3),
                   "Stillbirth"=sum(trans==4)),
             by=list(LM)]
nevent.m <- melt(nevent,id.vars = c("LM")) 

nriskevents<-rbindlist(list(nrisk.m,nevent.m),idcol="type")

plrisk<-ggplot(data=nrisk.m, aes(x = LM, y = value, fill = variable)) + 
  geom_bar(stat="identity")+
  xlab("Landmark time points (gestational week)")+
  ylab("Frequency")+
  scale_fill_brewer(palette="Paired")+
  guides(fill=guide_legend(title = "Exposure"))+
  ggtitle(paste("Exposure values"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plevent<-ggplot(data=nevent.m, aes(x = LM, y = value, fill = variable)) + 
  geom_bar(stat="identity", position = "stack")+
  labs(y="Frequency")+
  xlab("Landmark time points (gestational week)")+
  ggtitle(paste("Events at horizon"))+
  scale_fill_manual(values=clrs2,name="Event", labels=expo.lbls)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pl_compositionLMs<-arrangeGrob(plrisk,plevent,nrow=1)
plot(pl_compositionLMs)

ggsave(paste0(wdimages,"Unseld_Figure2.tiff"),pl_compositionLMs,width=10,height=4,units = "in",dpi=300)




