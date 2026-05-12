
## Supplemental figure (S5): TIE bias under different parameters
## Naive, HDMA & PCR

#clear environment
rm(list = ls())

# !diagnostics off

#set working directory
# dir_wd <- 
# setwd(dir_wd)

#add libraries
library("ggplot2")
# library("plotROC")
library("stringr")
library("cowplot")
library("parallel")
# library("pROC")
library("RColorBrewer")
library("EpiModel")
library("plyr")
library("readxl")
library("tidyr")
library("ggpubr")
library("reshape2")
library("magrittr")
library("ggforce")
library("psych")
library("dplyr")

# assign colors to all methods
methods_list<-c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","HIMA Compositional","MedTest","LDM-med-freq","LDM-med-omni3","SparseMCMM")
colors.methods<-c("#E377C2","#FF8000FF","#FBC02D","#33FF00FF","#F8AFA6","#E61A33FF","#4DAF4A","#1AB2FFFF",
                  "#664CFFFF","#A65628") 
names(colors.methods)<-methods_list
pie(x=seq(1,length(methods_list),1),col=colors.methods,labels=methods_list)

methods_list_fig <- methods_list[methods_list %in% c("PCR","HDMA","Naive")]
colors.methods_fig<-colors.methods[(names(colors.methods) %in% c("PCR","HDMA","Naive"))]

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
load(paste0(dir_wd,"main/res_total.RData"))
load(paste0(dir_wd,"data_creation/simOutcome.RData"))

res_total<-merge(x=res_total,y=simOutcome,by=c("outcome_index"),all.x=TRUE,all.y=TRUE)
                                     
res_total$manipulation_name <- case_when(
  res_total$nMicrobe == 50 ~ "Microbe_size",
  res_total$nMicrobe == 100 & res_total$TIE %in% c(0.05, 0.1, 0.2, 1) ~ "TIE_size"
)

df_panel<-data.frame("panel"=seq(1,2),"manipulation_name"=c("TIE_size", "Microbe_size")) 
df_panel$manipulation_var<-c("TIE_lab", "Microbe_lab") #nPathways
df_panel$grouping<-c("method","method")


# Make variables and variable labels
res_total$TIE_lab<-factor(res_total$TIE,
                          levels=c(0.05,0.1,0.2,1),
                          labels=c("0.05","0.1","0.2","1"))

res_total$Microbe_lab<-factor(res_total$TIE,
                              levels=c(0.05,0.1,0.2,1),
                              labels=c("0.05","0.1","0.2","1"))

res_total$method<-factor(res_total$method,
                         levels=c("naive","pcr","hdma"),
                         labels=methods_list_fig)

res_total<-merge(res_total,df_panel,by="manipulation_name",all.x=TRUE,all.y=TRUE)


df_axis<-NULL
for(i in 1:nrow(df_panel)){
  sub<-res_total[which(res_total$panel==df_panel$panel[i]),]
  manipulation_var_i<-df_panel$manipulation_var[i]
  grouping_i<-df_panel$grouping[i]
  sub[,manipulation_var_i]<-factor(sub[,manipulation_var_i])
  sub$method<-factor(sub$method)
  
  if(grouping_i=="method"){
  # Make x-axis labels for grouping by method and then by sample size
  nManVar<-length(unique(sub[,manipulation_var_i]))
  df_x_seq <- crossing(
    unique(sub$method),
    unique(sub[,manipulation_var_i])
  )
  colnames(df_x_seq)<-c("method","man_var")
  df_x_seq$x_axis_seq<-NA
  count=1
  bracket_min<-NULL
  bracket_max<-NULL
  for(j in unique(df_x_seq$method)){
    df_x_seq$x_axis_seq[which(df_x_seq$method==j)]<-seq(count,count+(nManVar-1),1)
    bracket_min<-c(bracket_min,count)
    bracket_max<-c(bracket_max,count+(nManVar-1))
    count=count+(nManVar+1)
  }
  }
  
  if(grouping_i=="manipulation"){
    # Make x-axis labels for grouping by sample size/manipulation variable and then method
    nMethod<-length(unique(sub$method))
    df_x_seq <- crossing(
      unique(sub[,manipulation_var_i]),
      unique(sub$method)
    )
    colnames(df_x_seq)<-c("man_var","method")
    df_x_seq$x_axis_seq<-NA
    count=1
    bracket_min<-NULL
    bracket_max<-NULL
    for(j in unique(df_x_seq$man_var)){
      df_x_seq$x_axis_seq[which(df_x_seq$man_var==j)]<-seq(count,count+(nMethod-1),1)
      bracket_min<-c(bracket_min,count)
      bracket_max<-c(bracket_max,count+(nMethod-1))
      count=count+(nMethod+1)
    }
  }
  df_x_seq$panel<-i
  df_axis<-rbind(df_axis,df_x_seq)
}
rm(df_x_seq,sub,bracket_max,bracket_min)
df_axis<-merge(x=df_axis,y=df_panel[,c("panel","manipulation_var")],by=c("panel"),all.x=TRUE)

res_total$manipulation_var_value<-NA
for(i in 1:nrow(res_total)){
  res_total$manipulation_var_value[i]<-as.character(res_total[i,res_total$manipulation_var[i]])
}
colnames(df_axis)[3]<-"manipulation_var_value"
res_total<-merge(x=res_total,y=df_axis,by=c("panel","method","manipulation_var","manipulation_var_value"),all.x=TRUE)
res_total$scaled_error<-(res_total$est_TIE-res_total$empTIE_bin)/res_total$empTIE_bin

res_total$manipulation_var_value<-factor(res_total$manipulation_var_value)


res_total$manipulation_name_lab<-factor(res_total$manipulation_name,
                                        levels=c("TIE_size","Microbe_size"),
                                        labels=c("TIE effect size",
                                                 "Number of microbes"))

res_total$my_lab_2<-factor(res_total$manipulation_name,
                           levels=c("TIE_size","Microbe_size"),
                           labels=c("A) 100 microbes, 10 mediators;\nTIE effect size (0.05, 0.1, 0.2, 1)",
                                    "B) 50 microbes, 10 mediators;\nTIE effect size (0.05, 0.1, 0.2, 1)"))

# save(res_total,file=paste0(dir_wd, "main/res_total_plot.RData"))


ggp<-ggplot(data=res_total,aes(x=x_axis_seq,y=scaled_error,group=x_axis_seq))+
  geom_boxplot(aes(color=method),outlier.shape = NA)+
  labs(x="Within methods, variation in the facet variable",y="Scaled Error in Estimated Total Indirect Effect")+
  facet_wrap_paginate(~my_lab_2,
                      ncol=2,
                      nrow=1,
                      scales="free")+
  scale_x_continuous(breaks=res_total$x_axis_seq,labels=NULL,
                     expand = expansion(add = c(0.5, 1.2)))+
  geom_hline(data=res_total,aes(yintercept=0), linetype="solid", color = "black", size=1)+
  scale_color_manual(name="Methods:",values=colors.methods_fig)+my_theme+
  coord_cartesian(ylim=c(-1.5,1.5))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="bottom", panel.grid.major.x = element_blank(),
        strip.text = element_text(family="Helvetica",face="bold",size=10),
        strip.background = element_rect(fill="lightblue", colour="black",size=1),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"))) #0.1

ggp


save_plot(paste0(dir_wd,"Figure/figureS5.png"),ggp,dpi=300,base_height=6.5,base_width=8)

