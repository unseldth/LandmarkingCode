###Set path for input and output
wd<-dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(wd,"/GlobalVariables.R"), local = T)
source(paste0(programs.path,"HelpFunctionsAnaGyrase.R"), local = T)


theme_set(theme_bw(base_size=16))

LMlong<-readRDS(paste0(data.path,"/LMlong.Rds"))


trans.labs<-event.lbls[-4]
names(trans.labs)<-1:3

pd2<-position_dodge(0.4)

###################################################################
## Figure 3 Comparison of probabilities from all models ###########
###################################################################


probsMstate<-readRDS(file=paste0(wdresults,"probsMstate.Rds"))
probsCrude<-readRDS(file=paste0(wdresults,"probsCrude.Rds"))
probsSuper<-readRDS(file=paste0(wdresults,"probsSuper.Rds"))
models.lbls<-c("Crude LM","Super LM","Mstate")
names(models.lbls)<-1:3

probslistUnadj<-list(probsCrude,probsSuper,probsMstate)
probs.estsUnadj<-rbindlist(probslistUnadj,idcol="Model")

unique(probs.estsUnadj$trans)
str(probs.estsUnadj$exposure)
probs.estsUnadj[,"exposure":=factor(exposure,levels=c("0","1"),labels=c("No prior exposure","Prior exposure"))]
probs.estsUnadj<-probs.estsUnadj[trans!=4]




# direct comparison of the exposure groups
pl_unadj_comparison_expsoure<-ggplot(probs.estsUnadj[trans!=4], aes(x = LM, y = CIF,color=exposure)) +
  geom_errorbar(aes(ymax =upper, ymin = lower), width=0, position = pd2) +
  geom_point(size=1.5, position = pd2)+
  ylab("Predicted cumulative incidence at horizon")+
  xlab("Landmark (weeks)")+
  theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  scale_color_hue("Exposure at landmark")+
  facet_grid(trans~Model, scales="free_y", labeller = labeller(trans=trans.labs,
                                                                          Model=models.lbls))+
  theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="lightgrey"), panel.border = element_rect(color = "black",fill=NA),
        text = element_text(size = 16))
pl_unadj_comparison_expsoure
ggsave(paste0(wdimages,"Unseld_Figure3.tiff"),pl_unadj_comparison_expsoure,width=8,height=8,unit="in",dpi=300)


#################################
# Figure S1: Adjusted for age ###
#################################
probsAge<-readRDS(file=paste0(wdresults,"probsSuperAdj.Rds"))


pl<-ggplot(probsAge[trans!=4], aes(x = LM, y = CIF, colour=exposure)) +
  geom_errorbar(aes(ymax =upper, ymin = lower), width=0, position = pd2)+
  geom_point(position = pd2)+
  theme_minimal()+
  ylab("Predicted cumulative incidences at horizon")+
  xlab("Landmark (weeks)")+
  scale_color_hue("Exposure at landmark")+
  theme(legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  facet_grid(trans~Age, scales="free_y", labeller = labeller(trans=trans.labs))+
  theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="lightgrey"), panel.border = element_rect(color = "black",fill=NA),
        legend.position="bottom", legend.box = "horizontal",legend.direction = "horizontal",
        text = element_text(size = 16)
        )
  
pl
ggsave(paste0(wdimages,"Simulation/S1_ProbsSuperAdj.tiff"),pl,width=8,height=8,unit="in",dpi=300)






