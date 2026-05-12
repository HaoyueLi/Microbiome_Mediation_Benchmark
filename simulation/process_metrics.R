## This script pre-processes metric information for the project

dir_wd <- paste0(here::here(),"/") #working directory
setwd(dir_wd)
dir_output <- paste0(dir_wd,"/sim_data/")

#add libraries
library("readxl")
library("plyr")
library("sendmailR")

#load simulation data
load(paste0(dir_output, "main/simOutcome_features.RData"))
load(paste0(dir_output, "main/simOutcome.RData"))

#load list of methods
dfMethods<-as.data.frame(read_excel(paste0(dir_wd,methods_list),na="NA"))
row.names(dfMethods)<-dfMethods$name


## define variables that are needed for script
total_variables<-c("TE","DE","TIE","run_time","runStatus","method","outcome_index","nPC")
if(fPerformHypothesisTesting==TRUE){
  total_variables<-c(total_variables,c("TIE_var","TIE_pvalue","TIE_pvalue_freq","TIE_pvalue_pa","TIE_pvalue_tran","TIE_pvalue_omni3"))
}
res_total<-NULL
res_component<-NULL


## Combine results for each method
for (j in 1:nrow(dfMethods)) {
  folder_path <- file.path(dir_output, "results_metrics", dfMethods$name[j])
  
  # List all RData files
  rdata_files <- list.files(
    path = folder_path,
    pattern = "result_\\d+\\.RData",  # Matches result_<outcome_index>.RData
    full.names = TRUE
  )
  
  results_metrics <- list()
  for (i in 1:length(rdata_files)) {
    tmp <- get(load(rdata_files[i]))
    results_metrics[[i]] <- tmp
  }
  save(results_metrics,file=paste0(dir_output,"results_metrics/",dfMethods$name[j],".RData"))
}

# check how many has error
# sum(sapply(results_metrics, function(x) !is.null(x$runStatus) && !(x$runStatus == "Success")))


## Extract data for res_total
for(j in 1:nrow(dfMethods)){
  load(paste0(dir_output, "results_metrics/",dfMethods$name[j],".RData"))
  
  # Extract total (non-feature specific) metrics for each method (ex. TE, DE, etc.)
  res1 <- NULL
  res2 <- NULL
  
  # LDM-med 
  if (dfMethods$name[j] == "ldm") {
    total_var <- total_variables[-10]
    
    df1_list <- lapply(results_metrics, function(x) {
      sapply(total_var, function(var) eval(parse(text = paste0("x$", var))))
    })
    
    df1 <- as.data.frame(do.call(rbind, df1_list))
    colnames(df1)<-total_var
    df1$TIE_pvalue <- NA
    res1<-df1
  }
  
  # other methods
  else { 
    total_var <- total_variables[1:(length(total_variables) - 4)]
    
    df2_list <- lapply(results_metrics, function(x) {
      unlist(x[total_var], use.names = FALSE)
    })
    
    df2 <- as.data.frame(do.call(rbind, df2_list))
    colnames(df2) <- total_var
    df2$TIE_pvalue_freq <- NA
    df2$TIE_pvalue_pa <- NA
    df2$TIE_pvalue_tran <- NA
    df2$TIE_pvalue_omni3 <- NA
    res2<-df2
  }
  res_total <- rbind(res_total,res1,res2)
}

## post-processing of res_total
res_total<-data.frame(res_total,stringsAsFactors = FALSE)
rownames(res_total)<-NULL
drop_rows<-which(res_total$outcome_index=="NULL")
if(length(drop_rows)>0){
  res_total<-res_total[-drop_rows,]
}
colnames(res_total)[1:3]<-paste0("est_",colnames(res_total)[1:3])
res_total$outcome_index<-as.numeric(res_total$outcome_index)

# remove those with error status
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
  res_total$TIE_pvalue<-as.numeric(unlist(res_total$TIE_pvalue))
  res_total$TIE_pvalue_freq<-as.numeric(res_total$TIE_pvalue_freq)
  res_total$TIE_pvalue_pa<-as.numeric(res_total$TIE_pvalue_pa)
  res_total$TIE_pvalue_tran<-as.numeric(res_total$TIE_pvalue_tran)
  res_total$TIE_pvalue_omni3<-as.numeric(res_total$TIE_pvalue_omni3)
  
  #Determine significance of the TIE (note hima does not provide an estimate of the significance of the TIE)
  res_total$TIE_sign<-NA
  res_total$TIE_sign_freq<-NA
  res_total$TIE_sign_pa<-NA
  res_total$TIE_sign_tran<-NA
  res_total$TIE_sign_omni3<-NA
  
  # for CCMM
  var_rows<-which(!(is.na(res_total$TIE_var)))
  res_total$TIE_pvalue[var_rows]<-2*pnorm(abs(res_total$est_TIE[var_rows]/sqrt(res_total$TIE_var[var_rows])), lower.tail=FALSE)
  
  pvals_rows<-which(!(is.na(res_total$TIE_pvalue)))
  res_total$TIE_sign[pvals_rows]<-ifelse(res_total$TIE_pvalue[pvals_rows]<0.05,TRUE,FALSE)
  
  # for LDM-med
  pvals_rows2<-which(!(is.na(res_total$TIE_pvalue_freq)))
  res_total$TIE_sign_freq[pvals_rows2]<-ifelse(res_total$TIE_pvalue_freq[pvals_rows2]<0.05,TRUE,FALSE)
  res_total$TIE_sign_pa[pvals_rows2]<-ifelse(res_total$TIE_pvalue_pa[pvals_rows2]<0.05,TRUE,FALSE)
  res_total$TIE_sign_tran[pvals_rows2]<-ifelse(res_total$TIE_pvalue_tran[pvals_rows2]<0.05,TRUE,FALSE)
  res_total$TIE_sign_omni3[pvals_rows2]<-ifelse(res_total$TIE_pvalue_omni3[pvals_rows2]<0.05,TRUE,FALSE)
  rm(pvals_rows,pvals_rows2)
  
  ## other methods already have TIE significance embedded in method functions
}


res_total_errors$message<-NULL
for(j in unique(res_total_errors$method)){
  load(paste0(dir_output, "results_metrics/",j,".RData"))
  index<-res_total_errors[which(res_total_errors$method==j),"outcome_index"]
  for(i in 1:length(index)){
    index_i<-index[i]
    matching_idx <- which(sapply(results_metrics, function(x) {
      !is.null(x$outcome_index) && x$outcome_index == index_i
    }))
    res_tmp <- results_metrics[[matching_idx]]
    res_total_errors$message[which(res_total_errors$method==j & res_total_errors$outcome_index==index_i)]<-
      paste(res_tmp$message, sep="", collapse="")
  }
}
save(res_total, file = paste0(dir_output, "main/res_total.RData"))
save(res_total_errors, file = paste0(dir_output, "main/res_total_errors.RData")) 


###############################################
## Extract component (feature-specific) metrics for each method (ex. CIE, a's, b's)
for(j in 1:nrow(dfMethods)){
  df=NULL
  
  # calculate which outcome_index for runs that did not error out
  viPassIndex<-res_total$outcome_index[which(res_total$method==dfMethods$name[j])]
  
  # for methods that return CIE information
  if(dfMethods$fComponentPathways[j]==TRUE & length(viPassIndex)>0){
    load(paste0(dir_output, "results_metrics/",dfMethods$name[j],".RData"))
    outcome_indices <- sapply(results_metrics, function(x) x$outcome_index)
    
    if(dfMethods$functionName[j]=="run_naive" ){
      for(i in 1:length(viPassIndex)){
        result_pos <- which(outcome_indices == viPassIndex[i])
        res_metric <- results_metrics[[result_pos]]
        df<-data.frame(outcome_index=res_metric$outcome_index,
                       method=res_metric$method,
                       feature=row.names(res_metric$feature_effects),
                       res_metric$feature_effects[,c("CIE","a","b")],
                       CIE_pvalue=res_metric$feature_effects[,c("CIE_pvalue")],
                       CIE_qvalue_freq=NA,
                       CIE_qvalue_pa=NA,
                       CIE_qvalue_tran=NA,
                       CIE_qvalue_omni3=NA,
                       CIE_lower_bound=NA,
                       CIE_upper_bound=NA,
                       CIE_pval_Bon=NA,
                       CIE_qval_FDR=NA,
                       row.names=NULL)
        res_component<-rbind(res_component,df)
      }
    }
    
    else if(dfMethods$functionName[j]=="run_hima_gaussian" | dfMethods$functionName[j]=="run_hima_compositional"){
      for(i in 1:length(viPassIndex)){
        result_pos <- which(outcome_indices == viPassIndex[i])
        res_metric <- results_metrics[[result_pos]]
        df<-data.frame(outcome_index=res_metric$outcome_index,
                       method=res_metric$method,
                       feature=row.names(res_metric$feature_effects),
                       res_metric$feature_effects[,c("CIE","a","b")],
                       CIE_pvalue=res_metric$feature_effects[,c("CIE_pvalue")],
                       CIE_qvalue_freq=NA,
                       CIE_qvalue_pa=NA,
                       CIE_qvalue_tran=NA,
                       CIE_qvalue_omni3=NA,
                       CIE_lower_bound=NA,
                       CIE_upper_bound=NA,
                       CIE_pval_Bon=NA,
                       CIE_qval_FDR=res_metric$feature_effects[,c("CIE_qvalue")],
                       row.names=NULL)
        if(sum(df$CIE==0)==nrow(df)){
          df$CIE_pvalue<-1
          df$CIE_pval_Bon<-1
          df$CIE_qval_FDR<-1
        }
        res_component<-rbind(res_component,df)
      }
    }
    
    else if(dfMethods$functionName[j]=="run_hdma"){
      for(i in 1:length(viPassIndex)){
        result_pos <- which(outcome_indices == viPassIndex[i])
        res_metric <- results_metrics[[result_pos]]
        df<-data.frame(outcome_index=res_metric$outcome_index,
                       method=res_metric$method,
                       feature=row.names(res_metric$feature_effects),
                       res_metric$feature_effects[,c("CIE","a","b")],
                       CIE_pvalue=res_metric$feature_effects[,c("CIE_pvalue")],
                       CIE_qvalue_freq=NA,
                       CIE_qvalue_pa=NA,
                       CIE_qvalue_tran=NA,
                       CIE_qvalue_omni3=NA,
                       CIE_lower_bound=NA,
                       CIE_upper_bound=NA,
                       CIE_pval_Bon=res_metric$feature_effects[,c("CIE_pvalue_Bonferroni")],
                       CIE_qval_FDR=res_metric$feature_effects[,c("CIE_qvalue")],
                       row.names=NULL)
        if(sum(df$CIE==0)==nrow(df)){
          df$CIE_pvalue<-1
          df$CIE_pval_Bon<-1
          df$CIE_qval_FDR<-1
        }
        res_component<-rbind(res_component,df)
      }
    }
    
    else if(dfMethods$functionName[j]=="run_ccmm"){
      for(i in 1:length(viPassIndex)){
        result_pos <- which(outcome_indices == viPassIndex[i])
        res_metric <- results_metrics[[result_pos]]
        Z<-res_metric$feature_effects$CIE/sqrt(res_metric$feature_effects$CIE_var)
        df<-data.frame(outcome_index=res_metric$outcome_index,
                       method=res_metric$method,
                       feature=row.names(res_metric$feature_effects),
                       res_metric$feature_effects[,c("CIE","a","b")],
                       CIE_pvalue=(2*pnorm(abs(Z), lower.tail=FALSE)),
                       CIE_qvalue_freq=NA,
                       CIE_qvalue_pa=NA,
                       CIE_qvalue_tran=NA,
                       CIE_qvalue_omni3=NA,
                       CIE_lower_bound=NA,
                       CIE_upper_bound=NA,
                       CIE_pval_Bon=NA,
                       CIE_qval_FDR=NA,
                       row.names=NULL)
        res_component<-rbind(res_component,df)
      }
    } 
    
    else if(dfMethods$functionName[j]=="run_ldm") {
      for(i in 1:length(viPassIndex)){
        result_pos <- which(outcome_indices == viPassIndex[i])
        res_metric <- results_metrics[[result_pos]]
        df<-data.frame(outcome_index=res_metric$outcome_index,
                       method=res_metric$method,
                       feature=row.names(res_metric$feature_effects),
                       res_metric$feature_effects[,c("CIE","a","b")],
                       CIE_pvalue=NA,
                       CIE_qvalue_freq=res_metric$feature_effects[,"CIE_qvalue_freq"],
                       CIE_qvalue_pa=res_metric$feature_effects[,"CIE_qvalue_pa"],
                       CIE_qvalue_tran=res_metric$feature_effects[,"CIE_qvalue_tran"],
                       CIE_qvalue_omni3=res_metric$feature_effects[,"CIE_qvalue_omni3"],
                       CIE_lower_bound=NA,
                       CIE_upper_bound=NA,
                       CIE_pval_Bon=NA,
                       CIE_qval_FDR=NA,
                       row.names=NULL)
        res_component<-rbind(res_component,df)
      }
    }
    else{
      print("Process metrics code for method is missing")
    } 
  }
}


## adjust p-value
# for ccmm and naive
if(fPerformHypothesisTesting==TRUE){ 
  
  for(j in c("naive","ccmm")) {
    
    for(k in unique(res_component$outcome_index[which(res_component$method==j)])){
      vals<-res_component$CIE_pvalue[which(res_component$method==j & res_component$outcome_index==k)]
      res_component$CIE_pval_Bon[which(res_component$method==j & res_component$outcome_index==k)]<-p.adjust(vals, method = "bonferroni", n = length(vals))
      res_component$CIE_qval_FDR[which(res_component$method==j & res_component$outcome_index==k)]<- p.adjust(vals, method = "fdr", n = length(vals))
    }
  }
  
}


## post-processing of res_component
colnames(res_component)[4:6]<-paste0("est_",colnames(res_component)[4:6])

#prep feature information
feat<-NULL

for(i in 1:length(simOutcome_features)){
  
  hold<-data.frame(outcome_index=i,feature=simOutcome_features[[i]]$feat_char$feature,fSpikedFeature=TRUE,
                   percent_nonzero=simOutcome_features[[i]]$feat_char$percent_nonzero,
                   path_vals=simOutcome_features[[i]]$feat_char$path_vals,
                   strength_a=simOutcome_features[[i]]$feat_char$strength_a)
  feat<-rbind(feat,hold)
}

# merge in feature information
res_component<-merge(x=res_component,y=feat,by=c("outcome_index","feature"),all.x=TRUE)
res_component$fSpikedFeature<-ifelse(!is.na(res_component$fSpikedFeature),TRUE,FALSE)
save(res_component, file = paste0(dir_output, "main/res_component.RData"))



