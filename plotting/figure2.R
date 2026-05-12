
## Figure 2: TIE TPR & FPR

#clear environment
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

# assign colors to all methods
methods_list<-c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","HIMA Compositional","MedTest",
                "LDM-med-freq","LDM-med-omni3","SparseMCMM")
colors.methods<-c("#E377C2","#FF8000FF","#FBC02D","#33FF00FF","#F8AFA6","#E61A33FF","#4DAF4A","#1AB2FFFF",
                  "#664CFFFF","#A65628")

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

# ######################################################################################################################
### PLOT 2 - TIE power and type I error

# Load simulation data
load(paste0(dir_output, "main/res_total.RData"))
load(paste0(dir_output, "data_creation/simOutcome.RData"))

## Need to reorganize data to have different LDM methods
## Make all methods only have one TIE p-value column
# reorganize other methods data
tmp1 <- res_total %>% dplyr::filter(!(method %in% c("ldm"))) %>%
  dplyr::select(est_TE,est_DE,est_TIE,run_time,runStatus,method,outcome_index,nPC,TIE_var,TIE_pvalue,TIE_sign) 

# reorganize ldm data
tmp2 <- res_total %>% dplyr::filter(method=="ldm") %>%
  dplyr::select(est_TE,est_DE,est_TIE,run_time,runStatus,method,outcome_index,nPC,TIE_var,TIE_pvalue_freq,TIE_pvalue_pa,TIE_pvalue_tran,TIE_pvalue_omni3,
                TIE_sign_freq,TIE_sign_pa,TIE_sign_tran,TIE_sign_omni3)
tmp2 <- tmp2 %>% tidyr::pivot_longer(
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

res_total_all <- rbind(tmp1, tmp2)
save(res_total_all,file = paste0(dir_output, "main/res_total_all.RData"))


simOutcome$TIE_grp<-NULL
for(i in 1:max(simOutcome$dataset_index)){
  rows<-which(simOutcome$dataset_index==i)
  simOutcome$TIE_grp[rows]<-max(simOutcome$TIE[rows])
}


## merge in true effect sizes
res_total_all<-merge(x=res_total_all,y=simOutcome,by="outcome_index",all.x=TRUE)
res_total_all<-res_total_all[-which(res_total_all$method=="hima_gaussian"),] #hima doesn't have an overall hypothesis test
res_total_all<-res_total_all[-which(res_total_all$method=="hima_compositional"),]
res_total_all<-res_total_all[-which(res_total_all$method=="hdma"),] #hdma doesn't have an overall hypothesis test
res_total_all<-res_total_all[-which(res_total_all$method=="ldm-pa"),]
res_total_all<-res_total_all[-which(res_total_all$method=="ldm-tran"),]

# customize the methods included in the figure
methods_list_TIE<-methods_list[!methods_list %in% c("HDMA","HIMA Gaussian","HIMA Compositional","SparseMCMM")]
colors.methods_TIE<-colors.methods[!(names(colors.methods) %in% c("HDMA","HIMA Gaussian","HIMA Compositional","SparseMCMM"))]


# Make facet labels
res_total_all$TIE_grp_lab<-factor(res_total_all$TIE_grp,
                                  levels=str_sort(unique(res_total_all$TIE_grp), numeric = TRUE,decreasing = FALSE),
                                  labels=c("Small TIE (0.5)","Medium TIE (1)", "Large TIE (1.5)"))

res_total_all$effect_type_lab<-factor(res_total_all$effect_type,
                                      levels=c("spiked","null_a","null_b","null_ab"),
                                      labels=c("Spiked","Null_a","Null_b","Null_ab"))

res_total_all$method<-factor(res_total_all$method,
                             levels=c("naive","pcr","ccmm","medtest","ldm-freq","ldm-omni3"),
                             labels=methods_list_TIE)


# Make x-axis labels for grouping by method and then by sample size
nSampleSizes<-length(unique(res_total_all$nSample))
df_x_seq <- crossing(
  unique(res_total_all$method),
  unique(res_total_all$nSample)
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
res_total_all<-merge(x=res_total_all,y=df_x_seq,by=c("method","nSample"),all.x=TRUE)


# MAKE TIE sensitivity, speficity plot
aggdata<-aggregate(TIE_sign ~ method+nSample+effect_type_lab+TIE_grp_lab, data=res_total_all, FUN=sum)
num<-aggregate(TIE_sign ~ method+nSample+effect_type_lab+TIE_grp_lab, data=res_total_all, FUN=length)
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
                    "method"=c("CCMM","CCMM","CCMM","CCMM"),stringsAsFactors=FALSE,
                    "nSample_lab"=c(1,1,1,1))

dummy$effect_type_lab_new<-factor(dummy$effect_type_lab_new,
                                  levels=levels(prop$effect_type_lab_new))
dummy$TIE_grp_lab<-factor(dummy$TIE_grp_lab,
                          levels=levels(prop$TIE_grp_lab))


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
  scale_color_manual(name="Methods:",values=colors.methods_TIE,guide = guide_legend())+my_theme+
  scale_x_continuous(breaks=df_x_seq$x_axis_seq,labels=df_x_seq$nSample)+
  # coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_text(family="Helvetica",size=6,angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major.x = element_blank())+
  theme(legend.position="bottom")
ggp

save_plot(paste0(dir_output,"Figure/TIE_figure/TIE_power.png"),ggp,dpi=300,base_height=6.5,base_width=(6.5))
ggsave(filename = paste0(dir_output,"Figure/TIE_figure/TIE_power.jpeg"),
       plot = ggp,
       dpi = 300,
       width = 6.5,
       height = 6.5,
       units = "in")

# save(prop,file = paste0(dir_output, "main/prop.RData"))
