
## Supplemental figures (S1&S2) of LDM-med variants

rm(list = ls())

# set working directory
# dir_wd <- 
# setwd(dir_wd)
# dir_output <- 


# #add libraries
library("ggplot2")
library("tidyr")
library("cowplot")
library("stringr")
library("reshape2")
library("parallel")
library("paletteer")

methods_list<-c("LDM-med-freq","LDM-med-omni3", "LDM-med-pa","LDM-med-tran")
colors.methods<-c("#1AB2FFFF","#664CFFFF","#984EA3", "#FF99BFFF" )

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
                           axis.line = element_line(colour = 'black', linewidth = 0.5),
                           axis.ticks = element_line(colour = "black", linewidth = 0.5),
                           strip.background = element_rect(fill="lightblue", colour="black",linewidth=1),
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

########################################################################################################################
load(paste0(dir_output, "main/res_total.RData"))
load(paste0(dir_output, "data_creation/simOutcome.RData"))

## Need to reorganize data to have different LDM methods
ldm_all <- res_total %>% dplyr::filter(method=="ldm") %>%
  dplyr::select(est_TE,est_DE,est_TIE,run_time,runStatus,method,outcome_index,nPC,TIE_var,TIE_pvalue_freq,TIE_pvalue_pa,TIE_pvalue_tran,TIE_pvalue_omni3,
                TIE_sign_freq,TIE_sign_pa,TIE_sign_tran,TIE_sign_omni3)
ldm_all <- ldm_all %>% tidyr::pivot_longer(
  cols = c(TIE_pvalue_freq, TIE_pvalue_pa, TIE_pvalue_tran, TIE_pvalue_omni3,
           TIE_sign_freq, TIE_sign_pa, TIE_sign_tran, TIE_sign_omni3),
  names_to = c(".value", "method_type"),
  names_pattern = "TIE_(pvalue|sign)_(freq|pa|tran|omni3)"
) %>%
  dplyr::mutate(
    method = paste0(method, "-", method_type)
  ) %>%
  dplyr::select(-method_type) %>%
  dplyr::rename(TIE_pvalue = pvalue, TIE_sign = sign)

# save(res_total_all,file = paste0(dir_output, "main/res_total_all.RData"))

simOutcome$TIE_grp<-NULL
for(i in 1:max(simOutcome$dataset_index)){
  rows<-which(simOutcome$dataset_index==i)
  simOutcome$TIE_grp[rows]<-max(simOutcome$TIE[rows])
}

# merge in true effect sizes
ldm_all<-merge(x=ldm_all,y=simOutcome,by="outcome_index",all.x=TRUE)


# Make facet labels
ldm_all$TIE_grp_lab<-factor(ldm_all$TIE_grp,
                                  levels=str_sort(unique(ldm_all$TIE_grp), numeric = TRUE,decreasing = FALSE),
                                  labels=c("Small TIE (0.5)","Medium TIE (1)", "Large TIE (1.5)"))

ldm_all$effect_type_lab<-factor(ldm_all$effect_type,
                                      levels=c("spiked","null_a","null_b","null_ab"),
                                      labels=c("Spiked","Null_a","Null_b","Null_ab"))

ldm_all$method<-factor(ldm_all$method,
                             levels=c("ldm-freq","ldm-omni3", "ldm-pa","ldm-tran"),
                             labels=methods_list)


# Make x-axis labels for grouping by method and then by sample size
nSampleSizes<-length(unique(ldm_all$nSample))
df_x_seq <- crossing(
  unique(ldm_all$method),
  unique(ldm_all$nSample)
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
ldm_all<-merge(x=ldm_all,y=df_x_seq,by=c("method","nSample"),all.x=TRUE)


aggdata<-aggregate(TIE_sign ~ method+nSample+effect_type_lab+TIE_grp_lab, data=ldm_all, FUN=sum)
num<-aggregate(TIE_sign ~ method+nSample+effect_type_lab+TIE_grp_lab, data=ldm_all, FUN=length)
colnames(aggdata)[5]<-"num_sign"
colnames(num)[5]<-"num_runs"
prop<-merge(x=num,y=aggdata,by=c("method","nSample","effect_type_lab","TIE_grp_lab"),all.x=TRUE,all.y=TRUE)
prop$p<-prop$num_sign/prop$num_runs
prop$se<-sqrt((prop$p*(1-prop$p))/prop$num_runs)


prop<-merge(x=prop,y=df_x_seq,by=c("method","nSample"),all.x=TRUE)
prop$effect_type_lab_new<-factor(prop$effect_type_lab,
                                 levels=c("Spiked","Null_a","Null_b","Null_ab"),
                                 labels=c("TPR","FPR (Null_a)","FPR (Null_b)", "FPR (Null_ab)"))


h_line_2<-tidyr::crossing("yintercept"=0.05,
                          "effect_type_lab_new"= c("FPR (Null_a)","FPR (Null_b)", "FPR (Null_ab)"),
                          "TIE_grp_lab"=levels(prop$TIE_grp_lab))
h_line_2$effect_type_lab_new<-factor(h_line_2$effect_type_lab_new,
                                     levels=c("FPR (Null_a)","FPR (Null_b)", "FPR (Null_ab)"))
h_line_2$TIE_grp_lab<-factor(h_line_2$TIE_grp_lab,
                             levels=levels(prop$TIE_grp_lab))


max(prop$p[which(prop$effect_type_lab=="Null_a")])
dummy <- data.frame(x_axis_seq = c(1,1,1,1), p = c(1,0.3,0.3,0.3),
                    "effect_type_lab_new" = levels(prop$effect_type_lab_new),
                    "TIE_grp_lab"=c("Small TIE (0.5)","Medium TIE (1)","Large TIE (1.5)","Small TIE (0.5)"),
                    "method"=c("LDM-med-freq","LDM-med-freq","LDM-med-freq","LDM-med-freq"),stringsAsFactors=FALSE,
                    "nSample_lab"=c(1,1,1,1))

dummy$effect_type_lab_new<-factor(dummy$effect_type_lab_new,
                                  levels=levels(prop$effect_type_lab_new))
dummy$TIE_grp_lab<-factor(dummy$TIE_grp_lab,
                          levels=levels(prop$TIE_grp_lab))


## figure S1 TIE
ggp<-ggplot(data=prop, aes(x=x_axis_seq, y=p, color=method)) +
  geom_line() +
  geom_point()+
  geom_hline(data=h_line_2,aes(yintercept = yintercept),linetype='dotted', linewidth=1)+
  labs(x="Sample Size",y="Value")+
  geom_blank(data=dummy) +
  facet_grid(cols=vars(TIE_grp_lab),rows=vars(effect_type_lab_new),
             scales="free_y",
             space="fixed",
             labeller="label_value")+
  scale_color_manual(name="Methods:",values=colors.methods,guide = guide_legend())+my_theme+
  scale_x_continuous(breaks=df_x_seq$x_axis_seq,labels=df_x_seq$nSample)+
  # coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(family="Helvetica",size=6,angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank())+
  theme(legend.position="bottom")
ggp

save_plot(paste0(dir_output,"Figure/TIE_figure/TIE_power_ldm.png"),ggp,dpi=300,
          base_height=6.5,base_width=(6.5))

ggsave(filename = paste0(dir_output,"Figure/TIE_figure/TIE_power_ldm.jpeg"),
       plot = ggp,
       dpi = 300,
       width = 6.5,
       height = 6.5,
       units = "in")

# save(prop,file = paste0(dir_output, "main/prop.RData"))


########################################################################################
## Figure S2 CIE
load(paste0(dir_output, "main/res_component.RData"))

ldm_component <- res_component %>% dplyr::filter(method=="ldm") %>%
  dplyr::select(-CIE_pvalue,-CIE_qval_FDR)
ldm_component <- ldm_component %>% tidyr::pivot_longer(
  cols = c(CIE_qvalue_freq, CIE_qvalue_pa, CIE_qvalue_tran, CIE_qvalue_omni3),
  names_to = c(".value", "method_type"),
  names_pattern = "(CIE_qvalue)_(freq|pa|tran|omni3)"
) %>%
  dplyr::mutate(
    method = paste0(method, "-", method_type)
  ) %>%
  dplyr::select(-method_type) %>%
  dplyr::rename(CIE_qval_FDR = CIE_qvalue)

# merge in true effect sizes
ldm_component_all<-merge(x=ldm_component,y=simOutcome,by="outcome_index",all.x=TRUE)


## Prepare to make plots of component indirect effects
ldm_component_all$CIE_bin<-ifelse(ldm_component_all$fSpikedFeature==TRUE & ldm_component_all$effect_type=="spiked",1,0)

# Make facet labels
ldm_component_all$TIE_grp_lab<-factor(ldm_component_all$TIE_grp,
                                      levels=str_sort(unique(ldm_component_all$TIE_grp), numeric = TRUE,decreasing = FALSE),
                                      labels=c("Small TIE (0.5)","Medium TIE (1)", "Large TIE (1.5)"))

ldm_component_all$effect_type_lab<-factor(ldm_component_all$effect_type,
                                          levels=c("spiked","null_a","null_b","null_ab"),
                                          labels=c("Spiked","Null_a","Null_b","Null_ab"))

ldm_component_all$method<-factor(ldm_component_all$method,
                                 levels=c("ldm-freq","ldm-omni3","ldm-pa","ldm-tran"),
                                 labels=methods_list)


# Make x-axis labels for grouping by method and then by sample size
nSampleSizes<-length(unique(ldm_component_all$nSample))
df_x_seq <- crossing(
  unique(ldm_component_all$method),
  unique(ldm_component_all$nSample)
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
ldm_component_all<-merge(x=ldm_component_all,y=df_x_seq,by=c("method","nSample"),all.x=TRUE)

ldm_component_all$CIE_qval_FDR[is.na(ldm_component_all$CIE_qval_FDR)]<-1


nonmed<-ldm_component_all %>%
  filter(CIE_bin==0) %>%
  group_by(method, outcome_index) %>%
  dplyr::summarise(num_nonmed=n(),
                   num_sign_nonmed_fdr_0.1 = sum(CIE_qval_FDR<0.1,na.rm=TRUE),
                   num_sign_nonmed_fdr_0.2 = sum(CIE_qval_FDR<0.2,na.rm=TRUE),
                   num_sign_nonmed_fdr_0.25 = sum(CIE_qval_FDR<0.25,na.rm=TRUE))

med<-ldm_component_all %>%
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
  dplyr::left_join(ldm_component_all %>% dplyr::select(outcome_index, TIE_grp_lab, effect_type_lab, x_axis_seq, method),
                   by=c("outcome_index", "method")) %>%
  dplyr::distinct()


my_list<-c("fdr_0.25") # can add more FDR thresholds

dummy <- data.frame("x_axis_seq" = c(1,1,1,1,1), FDR = c(0.4,0.5,0.5,0.5,0.5),
                    "effect_type_lab_new" = c("TPR","FDR (Spiked)","FDR (Null_a)","FDR (Null_b)","FDR (Null_ab)"),
                    "TIE_grp_lab"=c("Small TIE (0.5)","Medium TIE (1)","Large TIE (1.5)","Small TIE (0.5)","Small TIE (0.5)"),
                    "method"=c("LDM-med-freq","LDM-med-feq","LDM-med-freq","LDM-med-freq","LDM-med-freq"),stringsAsFactors=FALSE,
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
  
  # save(combs_long,file=paste0(dir_output,"main/combs_long.RData"))
  
  
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
    geom_hline(data=h_line_2,aes(yintercept = yintercept),linetype='dotted', linewidth=1)+
    geom_point() +
    geom_errorbar(aes(ymin = FDR - SE, ymax = FDR + SE)) +
    facet_grid(cols=vars(TIE_grp_lab),rows=vars(effect_type_lab_new),
               scales="free_y",
               space="fixed",
               labeller="label_value") +
    scale_color_manual(name="Methods:",values=colors.methods,guide = guide_legend())+
    my_theme+
    scale_x_continuous(breaks=df_x_seq$x_axis_seq,labels=df_x_seq$nSample)+
    theme(axis.text.x = element_text(family="Helvetica",size=6,angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major.x = element_blank())+
    theme(legend.position="bottom")+
    xlab("Sample Size")+ylab("Value")
  
  save_plot(paste0(dir_output,"Figure/CIE_figure/CIE_power_ldm.png"),plot,dpi=300,
            base_height=7.5,base_width=8.5)
  
  
}

