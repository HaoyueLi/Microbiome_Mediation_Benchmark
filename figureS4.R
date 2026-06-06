
## Figure S4: top CIE in MedDiet

library(dplyr)
library(ggplot2)

# set working directory
# dir_wd <- 
# setwd(dir_wd)

load(paste0(dir_wd,"results/res_component.RData"))
load(paste0(dir_wd,"results/res_component_sparsemcmm.RData"))

res_component_sparsemcmm$CIE_sign <- FALSE # no significant mediator at all
res_component_sparsemcmm$outcome_index <- ifelse(res_component_sparsemcmm$outcome_index==1, 
                                                 16, 17)

res_component <- rbind(res_component, res_component_sparsemcmm)

fdr_sign_level<-0.25

### check if any significant components in both datasets
## MedDiet
mlvs_results_main <- res_component %>% 
  filter(dataset=="MLVS" & outcome=="mscore_norm")
sum(mlvs_results_main$CIE_qval_FDR<fdr_sign_level, na.rm=TRUE) #1

sig_mlvs <- mlvs_results_main[mlvs_results_main$CIE_qval_FDR<fdr_sign_level 
                              & !is.na(mlvs_results_main$CIE_qval_FDR),]
# by HIMA Gaussian, Lactobacillus

mlvs_results_res <- res_component %>%
  filter(dataset=="MLVS" & outcome=="mscore_norm_residuals")
sum(mlvs_results_res$CIE_qval_FDR<fdr_sign_level, na.rm=TRUE) #0

## Blueberry RCT
blueberry_results <- res_component %>%
  filter(dataset=="blueberry")

sum(blueberry_results$CIE_qval_FDR < fdr_sign_level,na.rm=TRUE) #0
# no significant CIE detected in blueberry RCT


### Focus on MLVS dataset 
## filter features that have the largest CIE from at least 2 methods

load(paste0(dir_wd,"blueberry/names.RData"))
names_long$dataset<-"blueberry"
names<-names_long
load(paste0(dir_wd,"MLVS/names.RData"))
names_long$dataset<-"MLVS"
names<-rbind(names,names_long)
rm(names_long)

methods_list<-c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","HIMA Compositional","MedTest","LDM-med-freq","LDM-med-omni3","SparseMCMM")
colors.methods<-c("#E377C2","#FF8000FF","#FBC02D","#33FF00FF","#F8AFA6","#E61A33FF","#4DAF4A","#1AB2FFFF",
                  "#664CFFFF","#A65628") 
names(colors.methods)<-methods_list

methods_list_fig <- methods_list[methods_list %in% c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","SparseMCMM")]
colors.methods_fig<-colors.methods[(names(colors.methods) %in% c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","SparseMCMM"))]

shapes.methods_fig<-c(21,22,23,24,25,21)
names(shapes.methods_fig)<-methods_list_fig


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

#############################

names_genus<-unique(names[,-which(colnames(names) %in% c("species","long"))])

res_component<-res_component[which(res_component$outcome=="mscore_norm"|res_component$outcome=="cGMP_change"),]
res_component<-res_component[-which(res_component$method=="hima_compositional"),]
res_component<-merge(x=res_component,y=names_genus,by.x=c("dataset","feature"),by.y=c("dataset","genus"),all.x=TRUE)

res_component$method<-factor(res_component$method,
                              levels=c("naive","pcr","hdma","hima_gaussian","ccmm","SparseMCMM"),
                              labels=c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","SparseMCMM"))

res_component$my_color<-as.character(res_component$method)
res_component$my_color[which(res_component$CIE_sign)]<-"sign"

colors.methods.sign<-c(colors.methods,"#000000")
names(colors.methods.sign)<-c(methods_list,"sign")

res_component$feature_lab<-as.character(res_component$feature)
res_component$feature_lab<-sapply(strsplit(res_component$feature_lab,"_"), `[`, 1)
res_component$feature_lab<-substr(res_component$feature_lab,1,15)


# Make plot of most extreme (top 5, bottom 5) CIEs for MLVS - effect size, a's and b's 
res_component<-res_component[which(res_component$outcome=="mscore_norm" & res_component$microbiome_feature_level=="g"),]

df_tb<-NULL
top_x<-5
for(i in 1:length(unique(res_component$method))){
  res_component_i<-res_component[which(res_component$method==unique(res_component$method)[i]),]
  res_component_i<-res_component_i[order(res_component_i$est_CIE,decreasing=TRUE),]
  df_tb<-rbind(df_tb,(res_component_i[1:top_x,]))
  df_tb<-rbind(df_tb,(res_component_i[(nrow(res_component_i)-(top_x-1)):nrow(res_component_i),]))
}
df_tb$feature<-as.character(df_tb$feature)
tab<-table(df_tb$feature)
feat<-names(tab[which(tab>=2)])

res_component<-res_component[which(res_component$feature %in% feat),]

df_tb_long<-reshape::melt(res_component,id.vars=c("method","outcome_index","feature"),measure.vars=c("est_CIE","est_a","est_b"))
df_tb_long<-merge(x=df_tb_long,y=res_component[,c("method","feature","CIE_sign")],by=c("method","feature"),all.x=TRUE)
df_tb_long$my_color<-as.character(df_tb_long$method)
df_tb_long$my_color<-ifelse(df_tb_long$CIE_sign,"sign",df_tb_long$my_color)
df_tb_long<-merge(x=df_tb_long,y=names_genus[,c("genus","phylum")],by.x="feature",by.y="genus",all.x=TRUE)
df_tb_long$phylum<-factor(df_tb_long$phylum)

df_tb_long$phylum_lab<-as.character(df_tb_long$phylum)
my_rows<-which(df_tb_long$phylum_lab!="Firmicutes"&df_tb_long$phylum_lab!="Proteobacteria")
# shorten names of those phylum
df_tb_long$phylum_lab[my_rows]<-paste0(substr(df_tb_long$phylum_lab[my_rows],start=1,stop=5),".")
df_tb_long$phylum_lab<-factor(df_tb_long$phylum_lab)

y_labels = c("est_CIE" = "CIE",
             "est_a"= "a coefficient",
             "est_b"= "b coefficient")


############

## main plot
ggp <- ggplot(data=df_tb_long,aes(x=feature,
                                  y=value,
                                  ymin=value,
                                  ymax=value,
                                  color=my_color,
                                  fill=method,
                                  shape=method)) +
  geom_pointrange(fatten=2,
                  aes(),
                  position = position_jitterdodge(jitter.width = 0.15, #0.05
                                                  jitter.height = 0,
                                                  dodge.width = 0.4, #0
                                                  seed = 1),
                  stroke=1.5,
                  size=1.2) + 
  geom_hline(yintercept=0, 
             lty=2,
             lwd=0.5) +  # add a dotted line at x=1 after flip
  # coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("")+ylab("Estimated values")+my_theme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(rows=vars(variable),cols=vars(phylum_lab),scales="free",space="free_x",switch="y",
             labeller = labeller(variable = y_labels))+
  scale_fill_manual(name="Methods:",values=colors.methods_fig)+
  scale_shape_manual(name="Methods:",values=shapes.methods_fig)+
  scale_color_manual(name="Methods:",values=colors.methods.sign)+
  scale_y_continuous(position = "right")+
  theme() #legend.position = "bottom")
ggp

cowplot::save_plot(paste(dir_wd,"figS4_main.png"),ggp,dpi=150,base_height=6.5,base_width=10)


## for methods legend (shape+color)
ggp2 <- ggplot(data=df_tb_long,aes(x=feature,
                                  y=value,
                                  ymin=value,
                                  ymax=value,
                                  color=method,
                                  fill=method,
                                  shape=method)) +
  geom_pointrange(fatten=2,
                  aes(),
                  position = position_jitterdodge(jitter.width = 0.15, #0.05
                                                  jitter.height = 0,
                                                  dodge.width = 0.4, #0
                                                  seed = 1),
                  stroke=2) + 
  geom_hline(yintercept=0, 
             lty=2,
             lwd=0.5) +  # add a dotted line at x=1 after flip
  # coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("")+ylab("Estimated values")+my_theme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(rows=vars(variable),cols=vars(phylum_lab),scales="free",space="free_x",switch="y",
             labeller = labeller(variable = y_labels))+
  scale_fill_manual(name="Methods:",values=colors.methods_fig)+
  scale_shape_manual(name="Methods:",values=shapes.methods_fig)+
  scale_color_manual(name="Methods:",values=colors.methods_fig)+
  
  scale_y_continuous(position = "right")+
  theme()#legend.position = "bottom")
ggp2


cowplot::save_plot(paste(dir_wd,"figS4_legend.png"),ggp2,dpi=150,base_height=6.5,base_width=6.5)



## CIE significance legend
ggp3 <- ggp2 +  geom_point(
  data = subset(df_tb_long, CIE_sign == TRUE),
  aes(x = variable, y = value, shape = method, group = method, color = CIE_sign),
  fill = NA,
  stroke = 1,
  size = 3,
  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.05)
) +
  scale_color_manual(
    name = " ",
    values = c("TRUE" = "black"),
    labels = c("TRUE" = "CIE q<0.25"),
    na.translate = FALSE
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 21, fill = NA, stroke = 1.5, size = 4)
    )
  )

ggp3

cowplot::save_plot(paste(dir_wd,"figS4_legend_significance.png"),ggp3,dpi=150,base_height=6.5,base_width=6.5)

