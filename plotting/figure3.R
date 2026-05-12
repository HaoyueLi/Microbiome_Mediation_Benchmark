
## Figure 3: CIE TPR & FDR

#clear environment
rm(list = ls())

# set working directory
# dir_wd <- 
# setwd(dir_wd)
# dir_output <- 


# !diagnostics off

# #add libraries
library("ggplot2")
library("tidyr")
library("cowplot")
library("stringr")
library("reshape2")
library("dplyr")
library("paletteer")


my_theme<-theme_bw()+theme(axis.text.x = element_text(family="Helvetica",size=8),
                           axis.text.y = element_text(family="Helvetica",size=8),
                           legend.text=element_text(family="Helvetica",size=8),
                           legend.title=element_text(family="Helvetica",size=10),
                           axis.title.x = element_text(family="Helvetica",size=10),
                           axis.title.y = element_text(family="Helvetica",size=10),
                           strip.text = element_text(family="Helvetica",face="bold",size=10),
                           axis.line = element_line(colour = 'black', linewidth = 0.5),
                           axis.ticks = element_line(colour = "black", linewidth = 0.5),
                           strip.background = element_rect(fill="lightblue", colour="black",linewidth=1),
                           strip.text.x = element_text(margin = margin(.025, 0, .025, 0, "cm")),
                           panel.grid.minor = element_blank())

## assign colors to all methods
methods_list<-c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","MedTest","LDM-med-freq","LDM-med-omni3","SparseMCMM","HIMA Compositional")
colors.methods<-c("#E377C2","#FF8000FF","#FBC02D","#33FF00FF","#F8AFA6","#4DAF4A","#1AB2FFFF",
                  "#664CFFFF","#A65628","#E61A33FF") 
names(colors.methods)<-methods_list
pie(x=seq(1,length(methods_list),1),col=colors.methods,labels=methods_list)


# Load simulation data
load(paste0(dir_output, "main/res_total.RData"))
load(paste0(dir_output, "data_creation/simOutcome.RData"))
load(paste0(dir_output, "main/res_component.RData"))

## Need to reorganize data 
## Make all methods only have one TIE p-value column
tmp1 <- res_component %>% dplyr::filter(!(method %in% c("ldm"))) %>%
  dplyr::select(-CIE_pvalue,-CIE_qvalue_freq,-CIE_qvalue_pa,-CIE_qvalue_tran,-CIE_qvalue_omni3)

# ldm
tmp2 <- res_component %>% dplyr::filter(method=="ldm") %>%
  dplyr::select(-CIE_pvalue,-CIE_qval_FDR)
tmp2 <- tmp2 %>% tidyr::pivot_longer(
  cols = c(CIE_qvalue_freq, CIE_qvalue_pa, CIE_qvalue_tran, CIE_qvalue_omni3),
  names_to = c(".value", "method_type"),
  names_pattern = "(CIE_qvalue)_(freq|pa|tran|omni3)"
) %>%
  dplyr::mutate(
    method = paste0(method, "-", method_type)
  ) %>%
  dplyr::select(-method_type) %>%
  dplyr::rename(CIE_qval_FDR = CIE_qvalue)

res_component_all <- rbind(tmp1, tmp2)
save(res_component_all,file = paste0(dir_output, "main/res_component_all.RData"))


simOutcome$TIE_grp<-NULL
for(i in 1:max(simOutcome$dataset_index)){
  rows<-which(simOutcome$dataset_index==i)
  simOutcome$TIE_grp[rows]<-max(simOutcome$TIE[rows])
}

# merge in true effect sizes
res_component_all<-merge(x=res_component_all,y=simOutcome,by="outcome_index",all.x=TRUE)
res_component_all<-res_component_all[-which(res_component_all$method=="ldm-pa"),]
res_component_all<-res_component_all[-which(res_component_all$method=="ldm-tran"),]
res_component_all<-res_component_all[-which(res_component_all$method=="hima_compositional"),]


# customize the methods included in the figure
methods_list_CIE<-methods_list[!methods_list %in% c("MedTest","PCR","SparseMCMM", "HIMA Compositional")]
colors.methods_CIE<-colors.methods[!(names(colors.methods) %in% c("MedTest","PCR","SparseMCMM", "HIMA Compositional"))]

## Prepare to make plots of component indirect effects
res_component_all$CIE_bin<-ifelse(res_component_all$fSpikedFeature==TRUE & res_component_all$effect_type=="spiked",1,0)

# Make facet labels
res_component_all$TIE_grp_lab<-factor(res_component_all$TIE_grp,
                                      levels=str_sort(unique(res_component_all$TIE_grp), 
                                                      numeric = TRUE,decreasing = FALSE),
                                      labels=c("Small TIE (0.5)","Medium TIE (1)", 
                                               "Large TIE (1.5)"))

res_component_all$effect_type_lab<-factor(res_component_all$effect_type,
                                          levels=c("spiked","null_a","null_b","null_ab"),
                                          labels=c("Spiked","Null_a","Null_b","Null_ab"))

res_component_all$method<-factor(res_component_all$method,
                                 levels=c("naive","hdma","hima_gaussian","ccmm","ldm-freq","ldm-omni3"),
                                 labels=methods_list_CIE)


# Make x-axis labels for grouping by method and then by sample size
nSampleSizes<-length(unique(res_component_all$nSample))
df_x_seq <- crossing(
  unique(res_component_all$method),
  unique(res_component_all$nSample)
)
colnames(df_x_seq)<-c("method","nSample")
df_x_seq$x_axis_seq<-NA
count=1
bracket_min<-NULL
bracket_max<-NULL
for(i in unique(df_x_seq$method)){
  df_x_seq$x_axis_seq[which(df_x_seq$method==i)]<-seq(count,count+(nSampleSizes-1),1)
  bracket_min<-c(bracket_min,count)
  bracket_max<-c(bracket_max,count+(nSampleSizes-1))
  count=count+(nSampleSizes+1)
}
res_component_all<-merge(x=res_component_all,y=df_x_seq,by=c("method","nSample"),all.x=TRUE)

res_component_all$CIE_qval_FDR[is.na(res_component_all$CIE_qval_FDR)]<-1

# save(res_component_all,file=paste0(dir_output,"main/res_component_all_plot.RData"))


nonmed<-res_component_all %>%
  filter(CIE_bin==0) %>%
  group_by(method, outcome_index) %>%
  dplyr::summarise(num_nonmed=n(),
                   num_sign_nonmed_fdr_0.1 = sum(CIE_qval_FDR<0.1,na.rm=TRUE),
                   num_sign_nonmed_fdr_0.2 = sum(CIE_qval_FDR<0.2,na.rm=TRUE),
                   num_sign_nonmed_fdr_0.25 = sum(CIE_qval_FDR<0.25,na.rm=TRUE))

med<-res_component_all %>%
  filter(CIE_bin==1) %>%
  group_by(method, outcome_index) %>%
  dplyr::summarize(num_med=n(),
                   num_sign_med_fdr_0.1 = sum(CIE_qval_FDR<0.1,na.rm=TRUE),
                   num_sign_med_fdr_0.2 = sum(CIE_qval_FDR<0.2,na.rm=TRUE),
                   num_sign_med_fdr_0.25 = sum(CIE_qval_FDR<0.25,na.rm=TRUE))

combs<-merge(nonmed,med,by=c("outcome_index","method"),all.x=TRUE,all.y=TRUE)

## for all null scenarios, set the number of mediating features to 0
for(i in seq(3,10)){
  rows_na<-which(is.na(combs[,i]))
  combs[rows_na,i]<-0
}


## calculate power and FDR with q<0.25
combs$prop_sign_med_fdr_0.25<-combs$num_sign_med_fdr_0.25/(combs$num_med+0.00001)
combs$prop_sign_nonmed_fdr_0.25<-
  combs$num_sign_nonmed_fdr_0.25 / 
  (combs$num_sign_med_fdr_0.25 + combs$num_sign_nonmed_fdr_0.25+0.00001)

combs <- combs %>%
  dplyr::left_join(res_component_all %>% dplyr::select(outcome_index, TIE_grp_lab, effect_type_lab, x_axis_seq, method),
                   by=c("outcome_index", "method")) %>%
  dplyr::distinct()
# save(combs,file=paste0(dir_output,"main/combs.RData"))
# load(paste0(dir_output,"combs.RData"))


my_list<-c("fdr_0.25") # can add different FDR thresholds

dummy <- data.frame("x_axis_seq" = c(1,1,1,1,1), FDR = c(0.5,0.6,0.6,0.6,0.6),
                    "effect_type_lab_new" = c("TPR","FDR (Spiked)","FDR (Null_a)","FDR (Null_b)","FDR (Null_ab)"),
                    "TIE_grp_lab"=c("Small TIE (0.5)","Medium TIE (1)","Large TIE (1.5)","Small TIE (0.5)","Small TIE (0.5)"),
                    "method"=c("Naive","Naive","Naive","Naive","Naive"),stringsAsFactors=FALSE,
                    "nSample_lab"=c(1,1,1,1,1))


dummy$effect_type_lab_new<-factor(dummy$effect_type_lab_new,
                                  levels=c("TPR","FDR (Spiked)","FDR (Null_a)","FDR (Null_b)","FDR (Null_ab)"))
dummy$TIE_grp_lab<-factor(dummy$TIE_grp_lab,
                          levels=c("Small TIE (0.5)","Medium TIE (1)","Large TIE (1.5)"))



for(i in 1:length(my_list)){
  
  combs_long<-melt(combs[,c("outcome_index","method","TIE_grp_lab","effect_type_lab","x_axis_seq",
                            paste0("prop_sign_med_",my_list[i]),paste0("prop_sign_nonmed_",my_list[i]))], 
                   id=c("outcome_index","method","TIE_grp_lab","effect_type_lab","x_axis_seq"))
  combs_long<-combs_long[-which(combs_long$effect_type_lab!="Spiked" & combs_long$variable==paste0("prop_sign_med_",my_list[i])),]
  
  combs_long$effect_type_lab_new<-NULL
  combs_long$effect_type_lab_new[combs_long$effect_type_lab=="Spiked" & combs_long$variable==paste0("prop_sign_med_",my_list[i])]<-"TPR"
  combs_long$effect_type_lab_new[combs_long$effect_type_lab=="Spiked" & combs_long$variable==paste0("prop_sign_nonmed_",my_list[i])]<-"FDR (Spiked)"
  combs_long$effect_type_lab_new[combs_long$effect_type_lab=="Null_a"]<-"FDR (Null_a)"
  combs_long$effect_type_lab_new[combs_long$effect_type_lab=="Null_b"]<-"FDR (Null_b)"
  combs_long$effect_type_lab_new[combs_long$effect_type_lab=="Null_ab"]<-"FDR (Null_ab)"
  
  combs_long$effect_type_lab_new<-factor(combs_long$effect_type_lab_new,
                                         levels=c("TPR","FDR (Spiked)","FDR (Null_a)","FDR (Null_b)","FDR (Null_ab)"))
  

  h_line_2<-tidyr::crossing("yintercept"=0.25,
                            "effect_type_lab_new"= c("FDR (Spiked)","FDR (Null_a)","FDR (Null_b)", "FDR (Null_ab)"),
                            "TIE_grp_lab"=c("Small TIE (0.5)","Medium TIE (1)","Large TIE (1.5)"))
  h_line_2$effect_type_lab_new<-factor(h_line_2$effect_type_lab_new,
                                       levels=c("FDR (Spiked)","FDR (Null_a)","FDR (Null_b)", "FDR (Null_ab)"))
  h_line_2$TIE_grp_lab<-factor(h_line_2$TIE_grp_lab,
                               levels=c("Small TIE (0.5)","Medium TIE (1)","Large TIE (1.5)"))
 
  plot <- 
    combs_long %>% 
    dplyr::group_by(effect_type_lab_new, method, TIE_grp_lab, x_axis_seq) %>% 
    dplyr::summarise(FDR = mean(value),
                     SE = sd(value) / sqrt(dplyr::n())) %>%
    ggplot(aes(x = x_axis_seq, y = FDR, group=x_axis_seq,color=method)) +
    geom_blank(data=dummy) +
    geom_hline(data=h_line_2,aes(yintercept = yintercept),linetype='dotted', size=1)+
    geom_point() +
    geom_errorbar(aes(ymin = FDR - SE, ymax = FDR + SE)) +
    facet_grid(cols=vars(TIE_grp_lab),rows=vars(effect_type_lab_new),
               scales="free_y",
               space="fixed",
               labeller="label_value") +
    scale_color_manual(name="Methods:",values=colors.methods_CIE,guide = guide_legend())+
    my_theme+
    scale_x_continuous(breaks=df_x_seq$x_axis_seq,labels=df_x_seq$nSample)+
    theme(axis.text.x = element_text(family="Helvetica",size=6,angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major.x = element_blank())+
    theme(legend.position="bottom")+
    xlab("Sample Size")+ylab("Value")
  
  save_plot(paste0(dir_output,"/Figure/CIE_figure/CIE_power_",my_list[i],".png"),plot,dpi=300,base_height=7.5,base_width=8.5)
  
  ggsave(
    filename = paste0(dir_output,"/Figure/CIE_figure/CIE_power_",my_list[i],".pdf"),
    plot = plot,
    width = 8.5,
    height = 7.5,
    units = "in",
    dpi = 300
  )
  
}

