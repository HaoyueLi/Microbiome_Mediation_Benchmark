
## This script pre-processes metric information for mediation project
# except SparseMCMM and MedTest


#clear environment
rm(list = ls())

# set working directory
# dir_wd <- 
# setwd(dirwd)


methods_list<-"methods_vig.xlsx"
fPerformHypothesisTesting<-TRUE
fdr_sign_level<-0.25


#add libraries
library("readxl")
library("plyr")

dir_output<-paste0(dir_wd,"results/")

#load simulation data
load(paste0(dir_output, "res_vig.RData"))
load(paste0(dir_output, "df_vig.RData"))

#load list of methods
dfMethods<-as.data.frame(read_excel(paste0(dir_wd,methods_list),na="NA"))
row.names(dfMethods)<-dfMethods$name


#define variables that are needed for script
total_variables<-c("TE","DE","TIE","run_time","runStatus","method","outcome_index","nPC")
if(fPerformHypothesisTesting==TRUE){
  total_variables<-c(total_variables,c("TIE_var","TIE_pvalue"))
}
res_total<-NULL


## Extract total (non-feature specific) metrics for each method (ex. TE, DE, etc.)
df=NULL
for(i in 1:length(total_variables)){
  df<-cbind(df,sapply(res_vig, function(x) eval(parse(text=paste0("x$",total_variables[i])))))
}
res_total<-rbind(res_total,df)


## post-processing of res_total
res_total<-data.frame(res_total,stringsAsFactors = FALSE)
rownames(res_total)<-NULL
colnames(res_total)<-total_variables
colnames(res_total)[1:3]<-paste0("est_",colnames(res_total)[1:3])
res_total$outcome_index<-as.numeric(res_total$outcome_index)
res_total_errors<-NULL
if(length(which(res_total$runStatus=="Error"))>0){
  res_total_errors<-res_total[which(res_total$runStatus=="Error"),]
  res_total<-res_total[-which(res_total$runStatus=="Error"),]
}
res_total$est_TE<-as.numeric(res_total$est_TE)
res_total$est_DE<-as.numeric(res_total$est_DE)
res_total$est_TIE<-as.numeric(res_total$est_TIE)
res_total$run_time<-as.numeric(res_total$run_time)
res_total$nPC<-as.numeric(res_total$nPC)


if(fPerformHypothesisTesting==TRUE){
  res_total$TIE_var<-as.numeric(res_total$TIE_var)
  res_total$TIE_pvalue<-as.numeric(res_total$TIE_pvalue)

  #Determine significance of the TIE (note hima does not provide an estimate of the significance of the TIE)
  # for CCMM
  res_total$TIE_sign<-NA
  var_rows<-which(!(is.na(res_total$TIE_var)))
  res_total$TIE_pvalue[var_rows]<-2*pnorm(abs(res_total$est_TIE[var_rows]/sqrt(res_total$TIE_var[var_rows])), lower.tail=FALSE)
  pvals_rows<-which(!(is.na(res_total$TIE_pvalue)))
  res_total$TIE_sign[pvals_rows]<-ifelse(res_total$TIE_pvalue[pvals_rows]<0.05,TRUE,FALSE)
  rm(var_rows,pvals_rows)
}


######################

res_component<-NULL

## Extract component (feature-specific) metrics for each method (ex. CIE, a's, b's)
for(i in 1:length(res_vig)){
  df=NULL
  
  if(dfMethods[res_vig[[i]]$method,"fComponentPathways"]==TRUE & res_vig[[i]]$runStatus!="Error"){

    if(res_vig[[i]]$method=="naive" & fPerformHypothesisTesting==FALSE|res_vig[[i]]$method=="pcr"){
      df<-data.frame(outcome_index=res_vig[[i]]$outcome_index,
                     method=res_vig[[i]]$method,
                     feature=row.names(res_vig[[i]]$feature_effects),
                     res_vig[[i]]$feature_effects[,c("CIE","a","b")],
                     CIE_pvalue=res_vig[[i]]$feature_effects[,c("CIE_pvalue")],
                     CIE_lower_bound=NA,
                     CIE_upper_bound=NA,
                     CIE_pval_Bon=NA,
                     CIE_qval_FDR=NA,
                     row.names=NULL)
      res_component<-rbind(res_component,df)
    }
    
    else if(res_vig[[i]]$method=="naive" & fPerformHypothesisTesting==TRUE){
      df<-data.frame(outcome_index=res_vig[[i]]$outcome_index,
                     method=res_vig[[i]]$method,
                     feature=row.names(res_vig[[i]]$feature_effects),
                     res_vig[[i]]$feature_effects[,c("CIE","a","b")],
                     CIE_pvalue=res_vig[[i]]$feature_effects[,"CIE_pvalue"],
                     CIE_lower_bound=NA,
                     CIE_upper_bound=NA,
                     CIE_pval_Bon=NA,
                     CIE_qval_FDR=NA,
                     row.names=NULL)
      res_component<-rbind(res_component,df)
    }
    
    else if(res_vig[[i]]$method=="hima_gaussian"|res_vig[[i]]$method=="hima_compositional" ){
      df<-data.frame(outcome_index=res_vig[[i]]$outcome_index,
                     method=res_vig[[i]]$method,
                     feature=row.names(res_vig[[i]]$feature_effects),
                     res_vig[[i]]$feature_effects[,c("CIE","a","b")],
                     CIE_pvalue=res_vig[[i]]$feature_effects[,c("CIE_pvalue")],
                     CIE_lower_bound=NA,
                     CIE_upper_bound=NA,
                     CIE_pval_Bon=NA,
                     CIE_qval_FDR=res_vig[[i]]$feature_effects[,c("CIE_qvalue")],
                     row.names=NULL)
      if(sum(df$CIE==0)==nrow(df)){
        df$CIE_pval_Bon<-1
        df$CIE_qval_FDR<-1
      }
      res_component<-rbind(res_component,df)
    }
    
    else if(res_vig[[i]]$method=="hdma"){
      df<-data.frame(outcome_index=res_vig[[i]]$outcome_index,
                     method=res_vig[[i]]$method,
                     feature=row.names(res_vig[[i]]$feature_effects),
                     res_vig[[i]]$feature_effects[,c("CIE","a","b")],
                     CIE_pvalue=res_vig[[i]]$feature_effects[,c("CIE_pvalue")],
                     CIE_lower_bound=NA,
                     CIE_upper_bound=NA,
                     CIE_pval_Bon=NA,
                     CIE_qval_FDR=res_vig[[i]]$feature_effects[,c("CIE_qvalue")],
                     row.names=NULL)
      if(sum(df$CIE==0)==nrow(df)){
        df$CIE_pvalue<-1
        df$CIE_qval_FDR<-1
      }
      res_component<-rbind(res_component,df)
    }
    
    else if(res_vig[[i]]$method=="ccmm"){
      df<-data.frame(outcome_index=res_vig[[i]]$outcome_index,
                     method=res_vig[[i]]$method,
                     feature=row.names(res_vig[[i]]$feature_effects),
                     res_vig[[i]]$feature_effects[,c("CIE","a","b")],
                     CIE_pvalue=res_vig[[i]]$feature_effects[,c("CIE_pvalue")],
                     CIE_lower_bound=NA,
                     CIE_upper_bound=NA,
                     CIE_pval_Bon=NA,
                     CIE_qval_FDR=NA,
                     row.names=NULL)
      res_component<-rbind(res_component,df)
    } 
    
    else{
      print("Process metrics code for method is missing")
    } 
  }
}


## process significance testing data
# for CCMM and naive
if(fPerformHypothesisTesting==TRUE){ 
  
  for(j in c("naive","ccmm")) {
    for(k in unique(res_component$outcome_index[which(res_component$method==j)])){
      vals<-res_component$CIE_pvalue[which(res_component$method==j & res_component$outcome_index==k)]
      res_component$CIE_pval_Bon[which(res_component$method==j & res_component$outcome_index==k)]<-p.adjust(vals, method = "bonferroni", n = length(vals))
      res_component$CIE_qval_FDR[which(res_component$method==j & res_component$outcome_index==k)]<- p.adjust(vals, method = "fdr", n = length(vals))
    }
  }
  
  non_zero_rows<-which(!(is.na(res_component$CIE_qval_FDR)))
  res_component$CIE_sign<-NA
  res_component$CIE_sign[non_zero_rows]<-ifelse(res_component$CIE_qval_FDR[non_zero_rows]<fdr_sign_level,TRUE,FALSE)
}

## post-processing of res_component
colnames(res_component)[4:6]<-paste0("est_",colnames(res_component)[4:6])

res_component<-merge(x=res_component,y=df_vig,by.x=c("method","outcome_index"),by.y=c("method","index"),all.x=TRUE)
res_total<-merge(x=res_total,y=df_vig,by.x=c("method","outcome_index"),by.y=c("method","index"),all.x=TRUE)

save(res_component, file = paste0(dir_output, "res_component.RData"))
save(res_total, file = paste0(dir_output, "res_total.RData"))





