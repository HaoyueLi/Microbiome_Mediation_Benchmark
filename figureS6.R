
## Figure S6: spiked vs. empirical TIE


library("ggplot2")
library("cowplot")
library("stringr")
library("ggforce")

# set working directory
# dir_wd<-
# setwd(dir_wd)

load(paste0(wd,"sim_data/data_creation/simOutcome.RData"))
df<-simOutcome
df<-df[which(df$nMicrobe==100),]
df<-df[which(df$TIE==1),]
# there are replicates of nSample=100 but with different param_set_new
df_100 <- df[which(df$param_set_new==1),]
df <- df[-which(df$nSample==100),]
df <- rbind(df_100, df)
df$nSample_lab<-factor(df$nSample,
                       levels=str_sort(unique(df$nSample), numeric = TRUE),
                       labels=paste0("Samples:",str_sort(unique(df$nSample), numeric = TRUE)))
aggdata <-aggregate(df$empTIE_bin, by=list(df$nSample_lab),
                    FUN=mean, na.rm=TRUE)

ggp<-ggplot(data=df,aes(x=empTIE_bin))+
  geom_histogram(bins=50)+
  geom_vline(aes(xintercept = TIE),color="red")+
  facet_wrap_paginate(~nSample_lab,
                      ncol=1,
                      nrow=4,
                      scales="fixed")+
  # scale_y_continuous(limits = c(0, 20)) +
  xlab("Empirical TIE values")+
  ylab("Count")
ggp

save_plot(paste0(wd,"Figure/Figure_S6.png"),ggp,dpi=300,base_height=(6.5),base_width=(3.25))

