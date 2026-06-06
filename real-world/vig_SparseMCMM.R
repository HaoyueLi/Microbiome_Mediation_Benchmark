
## This script runs SparseMCMM on real-world data
# run each dataset separately


#clear environment
rm(list = ls())

#set working directory
# dir_wd<-
# setwd(dir_wd)

#load libraries
library("stringr")
library("dplyr")
library("ape")
library("readr")
library("SparseMCMM")

outcomes_mlvs<-c("mscore_norm","mscore_norm_residuals")
outcomes_blueberry <- c("cGMP_change_norm")

df_vig<-tidyr::crossing(method=c("SparseMCMM"), 
                        microbiome_feature_level=c("g"),
                        dataset=c("MLVS","blueberry"),
                        outcome=c(outcomes_mlvs,outcomes_blueberry))
df_vig<-df_vig[-which(df_vig$dataset=="MLVS" & df_vig$outcome %in% outcomes_blueberry),]
df_vig<-df_vig[-which(df_vig$dataset=="blueberry" & df_vig$outcome %in% outcomes_mlvs),]
df_vig$index<-seq(1,nrow(df_vig),1)
df_vig$exposure_var<-NA
df_vig$exposure_var[which(df_vig$dataset=="MLVS")]<-"treatment"
df_vig$exposure_var[which(df_vig$dataset=="blueberry")]<-"treatment_1"
df_vig$exposure_var[which(df_vig$outcome=="cGMP_change"|df_vig$outcome=="cGMP_change_norm")]<-"treatment_2"

time_index <- 1


###################
## run each dataset separately 
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)


## Blueberry RCT
param1<-df_vig[1,]
my_bugs1 <- get(load(paste0(dir_wd, param1$dataset, "/bugs_", param1$microbiome_feature_level, "_pc.RData")))
df1 <- get(load(paste0(dir_wd, param1$dataset, "/metadata.RData")))

rows1 <- which(!is.na(df1[, param1$outcome]))
Y1 <- df1[rows1, param1$outcome]
A1 <- df1[rows1, param1$exposure_var]
df_m1 <- my_bugs1[rows1, ]

row.names(df_m1)<-NULL
M1<-as.matrix(df_m1,dimnames=NULL)


run_time<-system.time(rslt.SparseMCMM1<-SparseMCMM(Treatment=A1,
                                                   otu.com=M1,
                                                   outcome=Y1,
                                                   num.per=100,
                                                   n.split=10,
                                                   ncores=50))[time_index]

CIE1=as.numeric(rslt.SparseMCMM1$`Compontent-wise ME`["Mean",])
feature_effects1<-data.frame("CIE"=CIE1,"a"=rep(NA),"b"=rep(NA),
                             "CIE_pvalue"=rep(NA),row.names=colnames(df_m1))

TE1<-as.numeric(unname(rslt.SparseMCMM1$`Esitmated Causal Effects`[1,3]))
DE1<-as.numeric(unname(rslt.SparseMCMM1$`Esitmated Causal Effects`[1,1]))
TIE1<-as.numeric(unname(rslt.SparseMCMM1$`Esitmated Causal Effects`[1,2]))
TIE_pvalue1<-as.numeric(rslt.SparseMCMM1$Test[["OME"]])

res1<-list(TE=TE1,
           DE=DE1,
           TIE=TIE1,
           run_time=run_time,
           feature_effects=feature_effects1,
           TIE_var=NA,
           TIE_pvalue=TIE_pvalue1,
           nPC=NA,
           raw_result=rslt.SparseMCMM1,
           runStatus="Success",
           method="SparseMCMM",
           outcome_index=1,
           message=NA)
save(res1, file = paste0(dir_wd, "results/res_sparsemcmm_blueberry.RData"))


## MedDiet (MLVS)
param2 <- df_vig[2,]
my_bugs2 <- get(load(paste0(dir_wd, param2$dataset, "/bugs_", param2$microbiome_feature_level, "_pc.RData")))
df2 <- get(load(paste0(dir_wd, param2$dataset, "/metadata.RData")))

rows2 <- which(!is.na(df2[, param2$outcome]))
Y2 <- df2[rows2, param2$outcome]
A2 <- df2[rows2, param2$exposure_var]
df_m2 <- my_bugs2[rows2, ]


row.names(df_m2)<-NULL
M2<-as.matrix(df_m2,dimnames=NULL)


run_time2<-system.time(rslt.SparseMCMM2<-SparseMCMM(Treatment=A2,
                                                   otu.com=M2,
                                                   outcome=Y2,
                                                   num.per=100,
                                                   n.split=10,
                                                   ncores=50))[time_index]

CIE2=as.numeric(rslt.SparseMCMM2$`Compontent-wise ME`["Mean",])
feature_effects2<-data.frame("CIE"=CIE2,"a"=rep(NA),"b"=rep(NA),
                             "CIE_pvalue"=rep(NA),row.names=colnames(df_m2))

TE2<-as.numeric(unname(rslt.SparseMCMM2$`Esitmated Causal Effects`[1,3]))
DE2<-as.numeric(unname(rslt.SparseMCMM2$`Esitmated Causal Effects`[1,1]))
TIE2<-as.numeric(unname(rslt.SparseMCMM2$`Esitmated Causal Effects`[1,2]))
TIE_pvalue2<-as.numeric(rslt.SparseMCMM2$Test[["OME"]])

res2<-list(TE=TE2,
           DE=DE2,
           TIE=TIE2,
           run_time=run_time2,
           feature_effects=feature_effects2,
           TIE_var=NA,
           TIE_pvalue=TIE_pvalue2,
           nPC=NA,
           raw_result=rslt.SparseMCMM2,
           runStatus="Success",
           method="SparseMCMM",
           outcome_index=2,
           message=NA)
save(res2, file = paste0(dir_wd, "results/res_vig_sparsemcmm_meddiet.RData"))


## MedDiet residual outcome
# param3 <- df_vig[3,]
# my_bugs3 <- get(load(paste0(dir_wd, param3$dataset, "/bugs_", param3$microbiome_feature_level, "_pc.RData")))
# df3 <- get(load(paste0(dir_wd, param3$dataset, "/metadata.RData")))
# 
# rows3 <- which(!is.na(df3[, param3$outcome]))
# Y3 <- df3[rows3, param3$outcome]
# A3 <- df3[rows3, param3$exposure_var]
# df_m3 <- my_bugs3[rows3, ]
# 
# row.names(df_m3)<-NULL
# M3<-as.matrix(df_m3,dimnames=NULL)
# colnames(M3)<-paste0("V",seq(1,ncol(M3)))
# 
# run_time3<-system.time(rslt.SparseMCMM3<-SparseMCMM(Treatment=A3,
#                                                     otu.com=M3,
#                                                     outcome=Y3,
#                                                     num.per=100,
#                                                     n.split=10,
#                                                     ncores=50))[time_index]
# 
# CIE3=as.numeric(rslt.SparseMCMM3$`Compontent-wise ME`["Mean",])
# feature_effects3<-data.frame("CIE"=CIE3,"a"=rep(NA),"b"=rep(NA),
#                              "CIE_pvalue"=rep(NA),row.names=colnames(df_m3))
# 
# TE3<-as.numeric(unname(rslt.SparseMCMM3$`Esitmated Causal Effects`[1,3]))
# DE3<-as.numeric(unname(rslt.SparseMCMM3$`Esitmated Causal Effects`[1,1]))
# TIE3<-as.numeric(unname(rslt.SparseMCMM3$`Esitmated Causal Effects`[1,2]))
# TIE_pvalue3<-as.numeric(rslt.SparseMCMM3$result$Test[["OME"]])
# 
# res3<-list(TE=TE3,
#            DE=DE3,
#            TIE=TIE3,
#            run_time=run_time3,
#            feature_effects=feature_effects3,
#            TIE_var=NA,
#            TIE_pvalue=TIE_pvalue3,
#            nPC=NA,
#            raw_result=rslt.SparseMCMM3,
#            runStatus="Success",
#            method="SparseMCMM",
#            outcome_index=3,
#            message=NA)







