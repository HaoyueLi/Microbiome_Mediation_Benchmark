
## Supplemental figure S3: CCMM TIE Estimation

rm(list = ls())
# set working directory
# dir_wd <- 
# setwd(dir_wd)

library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("tidyr")
library("dplyr")

methods_list<-c("ccmm")
colors.methods<-c("#F8AFA6")
names(colors.methods)<-methods_list
pie(x=seq(1,length(methods_list),1),col=colors.methods,labels=methods_list)

# Define a theme to use for all the plots
my_theme<-theme_bw()+theme(axis.text.x = element_text(family="Helvetica",size=8),
                           axis.text.y = element_text(family="Helvetica",size=8),
                           legend.text=element_text(family="Helvetica",size=8),
                           legend.title=element_text(family="Helvetica",size=10),
                           axis.title.x = element_text(family="Helvetica",size=10),
                           axis.title.y = element_text(family="Helvetica",size=10),
                           strip.text = element_text(family="Helvetica",face="bold",size=10),
                           axis.line = element_line(colour = 'black', size = 0.5),
                           axis.ticks = element_line(colour = "black", size = 0.5),
                           strip.background = element_rect(fill="lightblue", colour="black",size=1),
                           strip.text.x = element_text(margin = margin(.025, 0, .025, 0, "cm")),
                           panel.grid.minor = element_blank())

customCeiling <- function(x, Decimals=1) {
  x2<-x*10^Decimals
  ceiling(x2)/10^Decimals
}

customFloor <- function(x, Decimals=1) {
  x2<-x*10^Decimals
  floor(x2)/10^Decimals
}


######################################################################################################################
load(paste0(dir_wd,"/sim_data/main/res_total.RData"))
load(paste0(dir_wd,"/sim_data/data_creation/simOutcome.RData"))
# load(paste0(dir_wd,"/res_total_errors.RData"))

res_total<-merge(x=res_total,y=simOutcome,by=c("outcome_index"),all.x=TRUE,all.y=TRUE)

res_total$facet_lab<-paste0(res_total$nMicrobe, " features")
res_total$facet_lab<-factor(res_total$facet_lab,
                            str_sort(unique(res_total$facet_lab),numeric=TRUE))

res_total$scaled_error<-(res_total$est_TIE-res_total$empTIE_bin)/res_total$empTIE_bin
res_total$nSample_lab<-factor(res_total$nSample,
                              levels=c(50,100,200,400))

ggp<-ggplot(data=res_total,aes(x=nSample_lab,y=scaled_error,group=nSample_lab))+
  geom_boxplot(aes(color=method),outlier.shape = NA)+
  labs(x="Sample Size",y="Scaled Error in Estimated Total Indirect Effect")+
  ggforce::facet_wrap_paginate(~facet_lab,
                      nrow=1,
                      scales="free")+
  scale_x_discrete()+
  geom_hline(data=res_total,aes(yintercept=0), linetype="solid", color = "black", size=1)+
  scale_color_manual(guide="none",values=colors.methods)+my_theme+
  coord_cartesian(ylim=c(-2,2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="bottom", panel.grid.major.x = element_blank(),
        strip.text = element_text(family="Helvetica",face="bold",size=10),
        strip.background = element_rect(fill="lightblue", colour="black",size=1),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm")))

save_plot(paste0(dir_wd,"Figure/figureS3.png"),ggp,dpi=300,base_height=6.5,base_width=6.5)


