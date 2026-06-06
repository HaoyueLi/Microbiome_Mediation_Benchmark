
## Figure 4: TIE estimation 

#clear environment
rm(list = ls())

# !diagnostics off

# set working directory
# dir_wd <- 
# setwd(dir_wd)

#add libraries
library("ggplot2")
library("stringr")
library("cowplot")
library("parallel")
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

## assign colors to all methods
methods_list<-c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","HIMA Compositional","MedTest","LDM-med-freq","LDM-med-omni3","SparseMCMM")
colors.methods<-c("#E377C2","#FF8000FF","#FBC02D","#33FF00FF","#F8AFA6","#E61A33FF","#4DAF4A","#1AB2FFFF",
                  "#664CFFFF","#A65628") 
names(colors.methods)<-methods_list
pie(x=seq(1,length(methods_list),1),col=colors.methods,labels=methods_list)

methods_list_fig <- methods_list[!methods_list %in% c("MedTest","SparseMCMM","LDM-med-freq","LDM-med-omni3")]
colors.methods_fig<-colors.methods[!(names(colors.methods) %in% c("MedTest","SparseMCMM","LDM-med-freq","LDM-med-omni3"))]

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

##########################################################################################################

## load non-compositional data results
load(paste0(dir_wd,"noncomp/data_creation/simOutcome.RData"))
simOutcome$manipulation_name<-"non_comp"
simOutcome$noncomp<-TRUE
load(paste0(dir_wd,"noncomp/main/res_total.RData"))
non_comp<-merge(x=res_total,y=simOutcome,by="outcome_index",all.x=TRUE,all.y=TRUE)
rm(simOutcome,res_total)


## load other results
load(paste0(dir_wd,"main/res_total.RData"))
load(paste0(dir_wd,"data_creation/simOutcome.RData"))

res_total<-merge(x=res_total,y=simOutcome,by=c("outcome_index"),all.x=TRUE,all.y=TRUE)
res_total$manipulation_name<-case_when(res_total$nSample %in% c(50,200,400) ~ "sample_size",
                                       res_total$TIE %in% c(0.5,1.5) ~"TIE_size",
                                       res_total$nPathways %in% c(5,20) ~"num_med",
                                       res_total$DE %in% c(1,3) ~"percent_med",
                                       res_total$dirPathways_a==-1|res_total$dirPathways_b==-1 ~"dir_med")

df_panel<-data.frame("panel"=seq(1:6),
                     "manipulation_name"=c("TIE_size","num_med","dir_med","sample_size","percent_med","non_comp"))
df_panel$manipulation_var<-c("TIE_lab","nPathways","dirPathways_lab","nSample","DE_lab","noncomp")
df_panel$grouping<-c("method","method","method","method","method","method")

res_total_all<-bind_rows(res_total,non_comp)
rm(non_comp,res_total,simOutcome)


rows_baseline<-which(res_total_all$nSample==100 & res_total_all$nPathways==10 & 
                       res_total_all$dirPathways_a==1 & res_total_all$dirPathways_b==1 & 
                       res_total_all$TIE==1 & res_total_all$DE==0 & is.na(res_total_all$manipulation_name))
baseline<-res_total_all[rows_baseline,]
res_total_all<-res_total_all[-rows_baseline,]
for(i in c("non_comp","sample_size","TIE_size","num_med","percent_med","dir_med")){
  sub<-baseline
  sub$manipulation_name<-i
  res_total_all<-bind_rows(res_total_all,sub)
}
rm(sub,baseline)
res_total_all$noncomp[which(is.na(res_total_all$noncomp))]<-FALSE

# Make variables and variable labels
res_total_all$dirPathways_a_lab<-ifelse(res_total_all$dirPathways_a==1,"+","-")
res_total_all$dirPathways_b_lab<-ifelse(res_total_all$dirPathways_b==1,"+","-")
res_total_all$dirPathways_lab<-factor(paste0(res_total_all$dirPathways_a_lab,res_total_all$dirPathways_b_lab),
                                  levels=c("++","--","+-","-+"),
                                  labels=c("++","--","+-","-+"))

res_total_all$TIE_lab<-factor(res_total_all$TIE,
                          levels=str_sort(unique(res_total_all$TIE), numeric = TRUE,decreasing = FALSE),
                          labels=c("-1","0.5","1", "1.5"))

res_total_all$percent_med<-res_total_all$TIE/(res_total_all$TE)

res_total_all$DE_lab<-NA
res_total_all$DE_lab<-factor(res_total_all$percent_med,
                             levels=c(0.25,0.5,1),
                             labels=c("25%","50%","100%"))

res_total_all$noncomp<-factor(res_total_all$noncomp,
                              levels=c(FALSE,TRUE),
                              labels=c("C","NC"))

res_total_all$nPathways<-factor(res_total_all$nPathways,
                              levels=str_sort(unique(res_total_all$nPathways), numeric = TRUE,decreasing = FALSE))

res_total_all$nSample<-factor(res_total_all$nSample,
                                levels=str_sort(unique(res_total_all$nSample), numeric = TRUE,decreasing = FALSE))

res_total_all$method<-factor(res_total_all$method,
                         levels=c("naive","pcr","hdma","hima_gaussian","ccmm","hima_compositional"),
                         labels=methods_list_fig)

res_total_all<-merge(res_total_all,df_panel,by="manipulation_name",all.x=TRUE,all.y=TRUE)


## remove CCMM for compositionality panel
## removce HIMA Compositional from all results
res_total_all<-res_total_all[-which(res_total_all$method=="CCMM" & res_total_all$manipulation_name=="non_comp"),] 
res_total_all<-res_total_all[-which(res_total_all$method=="HIMA Compositional"),]

df_axis<-NULL
for(i in 1:nrow(df_panel)){
  sub<-res_total_all[which(res_total_all$panel==df_panel$panel[i]),]
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

res_total_all$manipulation_var_value<-NA
for(i in 1:nrow(res_total_all)){
  res_total_all$manipulation_var_value[i]<-as.character(res_total_all[i,res_total_all$manipulation_var[i]])
}
colnames(df_axis)[3]<-"manipulation_var_value"
res_total_all<-merge(x=res_total_all,y=df_axis,by=c("panel","method","manipulation_var","manipulation_var_value"),all.x=TRUE)
res_total_all$scaled_error<-(res_total_all$est_TIE-res_total_all$empTIE_bin)/res_total_all$empTIE_bin

df_panel$ymin=NA
df_panel$ymax=NA
for(i in 1:nrow(df_panel)){
  quants<-quantile(res_total_all$scaled_error[which(res_total_all$panel==df_panel$panel[i])])
  df_panel$ymax[i]=quants[4]+1.5*(quants[4]-quants[2])
  df_panel$ymin[i]=quants[2]-1.5*(quants[4]-quants[2])
}

res_total_all$manipulation_var_value<-factor(res_total_all$manipulation_var_value)


res_total_all$manipulation_name_lab<-factor(res_total_all$manipulation_name,
                                            levels=c("sample_size","TIE_size","num_med",
                                                     "dir_med","percent_med","non_comp"),
                                            labels=c("Sample size", "TIE effect size", "Number of mediators",
                                                     "Direction of mediation","Percent of effect mediated","Data compositionality"))

res_total_all$my_lab_2<-factor(res_total_all$manipulation_name,
                               levels=c("non_comp","sample_size","TIE_size","num_med",
                                        "percent_med","dir_med"),
                               labels=c("A) Data compositionality\n(comp, non-comp)",
                                        "B) Sample size\n(50, 100, 200, 400)",
                                        "C) TIE effect size\n(0.5, 1, 1.5)",
                                        "D) Number of mediators\n(5, 10, 20)",
                                        "E) Percent of effect mediated\n(25%, 50%, 100%)",
                                        "F) Direction of mediation\n(++, --, +-, -+)"))

# save(res_total_all,file=paste0(dir_wd, "main/res_total_all_plot.RData"))


ggp<-ggplot(data=res_total_all,aes(x=x_axis_seq,y=scaled_error,group=x_axis_seq))+
  geom_boxplot(aes(color=method),outlier.shape = NA)+
  labs(x="Within methods, variation in the facet variable",y="Scaled Error in Estimated Total Indirect Effect")+
  facet_wrap_paginate(~my_lab_2,
                      ncol=3,
                      nrow=2,
                      scales="free")+
  scale_x_continuous(breaks=res_total_all$x_axis_seq,labels=NULL,
                     expand = expansion(add = c(0.5, 1.2)))+
  geom_hline(data=res_total_all,aes(yintercept=0), linetype="solid", color = "black", size=1)+
  scale_color_manual(name="Methods:",values=colors.methods_fig)+my_theme+
  coord_cartesian(ylim=c(-1.5,1.5))+
  # geom_text(data=ann_text, aes(label=my_lab,x = xval, y = yval,color = NULL,group= NULL),hjust = 0.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position="bottom", panel.grid.major.x = element_blank(),
        strip.text = element_text(family="Helvetica",face="bold",size=10),
        strip.background = element_rect(fill="lightblue", colour="black",size=1),
        strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"))) #0.1

ggp


save_plot(paste0(dir_wd,"Figure/fig4.png"),ggp,dpi=300,base_height=6.5,base_width=8)

