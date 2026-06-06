
## This script processes the results of SparseMCMM

#clear environment
rm(list = ls())

library("readxl")
library("plyr")

#set working directory
# dir_wd<-
# setwd(dir_wd)
dir_output<-paste0(dir_wd,"results/")

methods_list<-"methods_vig.xlsx"
fPerformHypothesisTesting<-TRUE
fdr_sign_level<-0.25


## metadata
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


## load results 
load(paste0(dir_output, "res_vig_sparsemcmm_blueberry.RData"))
load(paste0(dir_output, "res_vig_sparsemcmm_meddiet.RData"))

res_sparsemcmm <- list(res1, res2)

#define variables that are needed for script
total_variables<-c("TE","DE","TIE","run_time","runStatus","method","outcome_index","nPC")
if(fPerformHypothesisTesting==TRUE){
  total_variables<-c(total_variables,c("TIE_var","TIE_pvalue"))
}
res_total_sparsemcmm<-NULL


## Extract total (non-feature specific) metrics for each method (ex. TE, DE, etc.)
df=NULL
for(i in 1:length(total_variables)){
  df<-cbind(df,sapply(res_sparsemcmm, function(x) eval(parse(text=paste0("x$",total_variables[i])))))
}
res_total_sparsemcmm <- rbind(res_total_sparsemcmm, df)


## post-processing of res_total
res_total_sparsemcmm<-data.frame(res_total_sparsemcmm,stringsAsFactors = FALSE)
rownames(res_total_sparsemcmm)<-NULL
colnames(res_total_sparsemcmm)<-total_variables
colnames(res_total_sparsemcmm)[1:3]<-paste0("est_",colnames(res_total_sparsemcmm)[1:3])
res_total_sparsemcmm$outcome_index<-as.numeric(res_total_sparsemcmm$outcome_index)

res_total_sparsemcmm$est_TE<-as.numeric(res_total_sparsemcmm$est_TE)
res_total_sparsemcmm$est_DE<-as.numeric(res_total_sparsemcmm$est_DE)
res_total_sparsemcmm$est_TIE<-as.numeric(res_total_sparsemcmm$est_TIE)
res_total_sparsemcmm$run_time<-as.numeric(res_total_sparsemcmm$run_time)
res_total_sparsemcmm$nPC<-as.numeric(res_total_sparsemcmm$nPC)

res_total_sparsemcmm$TIE_sign<-ifelse(res_total_sparsemcmm$TIE_pvalue<0.05,TRUE,FALSE)

##########
res_component_sparsemcmm<-NULL

## Extract component (feature-specific) metrics for each method (ex. CIE, a's, b's)
for(i in 1:length(res_sparsemcmm)){
  df=NULL
  
  df <- data.frame(outcome_index=res_sparsemcmm[[i]]$outcome_index,
                   method=res_sparsemcmm[[i]]$method,
                   feature=row.names(res_sparsemcmm[[i]]$feature_effects),
                   res_sparsemcmm[[i]]$feature_effects[,c("CIE","a","b")],
                   CIE_pvalue=res_sparsemcmm[[i]]$feature_effects[,c("CIE_pvalue")],
                   CIE_lower_bound=NA,
                   CIE_upper_bound=NA,
                   CIE_pval_Bon=NA,
                   CIE_qval_FDR=NA,
                   row.names=NULL)
  res_component_sparsemcmm<-rbind(res_component_sparsemcmm,df)
}

colnames(res_component_sparsemcmm)[4:6]<-paste0("est_",colnames(res_component_sparsemcmm)[4:6])


res_component_sparsemcmm<-merge(x=res_component_sparsemcmm,y=df_vig,by.x=c("method","outcome_index"),by.y=c("method","index"),all.x=TRUE)
res_total_sparsemcmm<-merge(x=res_total_sparsemcmm,y=df_vig,by.x=c("method","outcome_index"),by.y=c("method","index"),all.x=TRUE)

save(res_component_sparsemcmm, file = paste0(dir_output, "res_component_sparsemcmm.RData"))
save(res_total_sparsemcmm, file = paste0(dir_output, "res_total_sparsemcmm.RData"))


