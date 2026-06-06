
## Figure 5: TE & TIE estimation 

rm(list = ls())

## set working directory
# dir_wd <- 
# setwd(dir_wd)

nCores=1



#load libraries
library("ggplot2")
library("cowplot")
library("vegan")
library("labdsv")
library("RColorBrewer")
library("RColorBrewer")
library("EpiModel")
library("tibble")
library("ggplot2")
library("viridis")
library("cowplot")
library("grid")
library("gridExtra")
library("dplyr")
library("Cairo")

#assign colors to all methods
methods_list<-c("Naive","PCR","HDMA","HIMA Gaussian","CCMM","HIMA Compositional","MedTest","LDM-med-freq","LDM-med-omni3","SparseMCMM")
colors.methods<-c("#E377C2","#FF8000FF","#FBC02D","#33FF00FF","#F8AFA6","#E61A33FF","#4DAF4A","#1AB2FFFF",
                  "#664CFFFF","#A65628") 
names(colors.methods)<-methods_list
pie(x=seq(1,length(methods_list),1),col=colors.methods,labels=methods_list)

methods_list_fig <- methods_list[methods_list %in% c("Naive","PCR","HDMA","HIMA Gaussian","CCMM", "SparseMCMM")]
colors.methods_fig<-colors.methods[(names(colors.methods) %in% c("Naive","PCR","HDMA","HIMA Gaussian","CCMM", "SparseMCMM"))]

shapes.methods_fig<-c(21,22,23,24,25,21) 
names(shapes.methods_fig)<-methods_list_fig


my_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(family = "Helvetica", size = 8),
    axis.text.y = element_text(family = "Helvetica", size = 8),
    legend.text = element_text(family = "Helvetica", size = 8),
    legend.title = element_text(family = "Helvetica", size = 10),
    axis.title.x = element_text(family = "Helvetica", size = 10),
    axis.title.y = element_text(family = "Helvetica", size = 10),
    strip.text = element_text(family = "Helvetica", face = "bold", size = 10),
    axis.ticks = element_line(colour = "black", size = 0.5),
    strip.background = element_rect(fill="lightblue", colour="black",linewidth=1),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank()  # <-- this removes extra overlap
  )

############################################################
#Load datasets
load(paste0(dir_wd,"results/res_total.RData"))
load(paste0(dir_wd,"results/res_total_sparsemcmm.RData"))

res_total_sparsemcmm$outcome_index <- c(19, 20)

res_total <- rbind(res_total, res_total_sparsemcmm)

# Analyze at the genus level for now
res_total_long<-reshape::melt(res_total,id.vars=c("method","outcome_index","dataset","outcome"),measure.vars=c("est_TIE","est_TE"))

res_total_long<-merge(x=res_total_long,y=res_total[,c("method","outcome_index","dataset","outcome","TIE_pvalue")],
                      by=c("method","outcome_index","dataset","outcome"),all.x=TRUE)
res_total_long$TIE_sign<-ifelse(res_total_long$TIE_pvalue<0.1 & res_total_long$variable=="est_TIE",TRUE,FALSE)

# exclude HIMA Compositional 
res_total_long<-res_total_long[-which(res_total_long$method=="hima_compositional"),]

res_total_long$method<-factor(res_total_long$method,
                              levels=c("naive","pcr","hdma","hima_gaussian","ccmm", "SparseMCMM"),
                              labels=c("Naive","PCR","HDMA","HIMA Gaussian","CCMM", "SparseMCMM"))

res_total_long$my_color<-as.character(res_total_long$method)
res_total_long$my_color[which(res_total_long$TIE_sign)]<-"sign"

colors.methods.sign<-c(colors.methods_fig,"#000000")
names(colors.methods.sign)<-c(methods_list_fig,"sign")


sub<-res_total_long[which(res_total_long$outcome=="mscore_norm"|res_total_long$outcome=="cGMP_change_norm"),]

# Label mapping for x-axis (later will be flipped to y axis)
x_labels <- c("est_TE" = "Estimated Total Effect", 
              "est_TIE" = "Estimated Total Indirect Effect")

sub$dataset_label <- factor(sub$dataset,
                            levels = c("MLVS", "blueberry"),
                            labels = c("MedDiet", "Blueberry RCT"))

sub$outcome_label <- factor(sub$outcome,
                            levels = c("mscore_norm", "cGMP_change_norm"),
                            labels = c("CMD biomarker score (normalized)",
                                       "Change in cGMP (normalized)"))


## main plot
ggp <- ggplot(data=sub,aes(x=variable,
                           y=value,
                           ymin=value,
                           ymax=value,
                           color=my_color,
                           fill=method,
                           shape=method,
                           group=method)) +
  geom_point(alpha=0.8,
             aes(),
             position = position_dodge(
               width = 0.6),
             stroke=1.5,
             size=3) +
  geom_hline(yintercept=0,
             lty=2,
             lwd=0.5) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("")+
  ylab("Estimated effect size")+
  my_theme+
  facet_wrap(~dataset_label+outcome_label,scales="free",
             nrow=1,
             ncol=2)+
  
  scale_fill_manual(name="Methods:",values=colors.methods_fig)+
  scale_shape_manual(name="Methods:",values=shapes.methods_fig)+
  scale_color_manual(name="Methods",values=colors.methods.sign)+
  
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous() 

ggp

save_plot(paste(dir_wd,"fig5_main.png"),ggp,dpi=400,base_height=4,base_width=10)
ggsave(filename = paste0(dir_wd, "fig5_main.jpeg"),
       plot = ggp,
       dpi = 400,
       width = 10,
       height = 4,
       units = "in")


## for methods legend (shape+color)
ggp2 <- ggplot(data=sub,aes(x=variable,
                            y=value,
                            ymin=value,
                            ymax=value,
                            color=method,
                            fill=method,
                            shape=method)) +
  geom_pointrange(fatten=5,stroke=1.5, #stroke=2 # fatten=5
                  aes(),
                  position = position_jitterdodge(jitter.width = 0.15,
                                                  jitter.height = 0,
                                                  dodge.width = 0.05,
                                                  seed = 1)) +
  geom_hline(yintercept=0,
             lty=2,
             lwd=0.5) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("")+
  ylab("Estimated effect size")+
  my_theme+
  facet_wrap(~dataset_label+outcome_label,scales="free",
             nrow=1,
             ncol=2)+
  scale_fill_manual(name="Methods:",values=colors.methods_fig)+
  scale_shape_manual(name="Methods:",values=shapes.methods_fig)+
  scale_color_manual(name="Methods",values=colors.methods_fig)+
  
  guides(
    color = guide_legend(
      override.aes = list(
        fill = colors.methods_fig,
        color = colors.methods_fig,
        stroke = 0.5,  # thin outline, or set to 0 to remove
        size=1
      )
    ),
    shape = guide_legend(
      override.aes = list(
        fill = colors.methods_fig,
        color = colors.methods_fig,
        stroke = 0.5,
        size=1
      )
    )
  ) +
  
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous() 

ggp2

save_plot(paste(dir_wd,"fig5_legend.png"),ggp2,dpi=400,base_height=4,base_width=10)
ggsave(filename = paste0(dir_wd, "fig5_legend.jpeg"),
       plot = ggp2,
       dpi = 400,
       width = 10,
       height = 4,
       units = "in")


## TIE significance legend
ggp3 <- ggp2 +  geom_point(
  data = subset(sub, TIE_sign == TRUE),
  aes(x = variable, y = value, shape = method, group = method, color = TIE_sign),
  fill = NA,
  stroke = 1,
  size = 3,
  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.05)
) +
  scale_color_manual(
    name = " ",
    values = c("TRUE" = "black"),
    labels = c("TRUE" = "Total indirect effect p<0.1"),
    na.translate = FALSE
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 21, fill = NA, stroke = 1.5, size = 3)
    )
  )

ggp3


save_plot(paste(dir_wd,"fig5_legend_sig.png"),ggp3,dpi=400,base_height=4,base_width=10)
ggsave(filename = paste0(dir_wd, "fig5_legend_sig.jpeg"),
       plot = ggp3,
       dpi = 400,
       width = 10,
       height = 4,
       units = "in")



## some stats calculation for results section
res_total_mlvs <- res_total[res_total$outcome=="mscore_norm" & res_total$method %in% c("ccmm","naive","pcr","SparseMCMM"),]
res_total_mlvs$est_TIE
res_total_mlvs$est_TE
res_total_mlvs$est_TIE/res_total_mlvs$est_TE

res_total_blueberry <- res_total[res_total$outcome=="cGMP_change_norm" & res_total$method %in% c("ccmm","naive","pcr","SparseMCMM") ,]
res_total_blueberry$est_TIE
res_total_blueberry$est_TE
res_total_blueberry$est_TIE/res_total_blueberry$est_TE
