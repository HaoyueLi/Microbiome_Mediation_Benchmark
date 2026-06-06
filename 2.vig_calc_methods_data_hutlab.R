
## This script contains methods on vignette datasets for the mediation project


#clear environment
rm(list = ls())

# #set working directory

# dir_wd<-
# setwd(dir_wd)

#load libraries
library("HIMA")
library("stringr")
library("ccmm")
# library("SparseMCMM")
library("readxl")
library("vegan")
source("helper_functions.R")
library("mvtnorm")
library("Matrix")
library("Rfast")
library("matrixStats")
library("tidyr")


#Run methods in example datasets
methods_list<-NA
fPerformHypothesisTesting<-TRUE
nIterationsHypothesisTesting<-500
fUseNonCompositionalData<-FALSE
time_index<-1 #which value of system.time to use

# To make sure parallel computing run using tmux, load df_vig here
load(paste0(dir_wd, "results/df_vig.RData"))

 
## Methods
run_ccmm<-function(A,Y,df_m,optionalParameters,outcome_index){
  M<-as.matrix(df_m)
  if(fPerformHypothesisTesting==FALSE){
    run_time<-system.time(rslt.ccmm <- ccmm(y=Y,M=M,tr=A,x=NULL))[time_index]
    feature_effects<-data.frame("CIE"=rslt.ccmm$IDEs,"a"=rep(NA),"b"=rep(NA),
                                "CIE_pvalue"=rep(NA),row.names=colnames(df_m))
    TIE_var<-NA
    TIE_pvalue<-NA
  }else{
    run_time<-system.time(rslt.ccmm <- ccmm(y=Y,M=M,tr=A,x=NULL,method.est.cov="normal"))[time_index]
    Z<-rslt.ccmm$IDEs/sqrt(rslt.ccmm$Var.IDEs)
    CIE_pvalue=(2*pnorm(abs(Z), lower.tail=FALSE))
    feature_effects<-data.frame("CIE"=rslt.ccmm$IDEs,"a"=rep(NA),"b"=rep(NA),
                                "CIE_pvalue"=CIE_pvalue,row.names=colnames(df_m))
    TIE_var<-unname(rslt.ccmm$Var.TIDE)
    TIE_pvalue<-2*pnorm(abs(rslt.ccmm$TIDE/sqrt(TIE_var)), lower.tail=FALSE)
  }
  DE<-rslt.ccmm$DE
  TIE<-rslt.ccmm$TIDE
  TE<-DE+TIE 
  return(list(TE=TE,
              DE=DE,
              TIE=TIE,
              run_time=run_time,
              feature_effects=feature_effects,
              TIE_var=TIE_var,
              TIE_pvalue=TIE_pvalue,
              nPC=NA))
}



run_hima_gaussian<-function(A,Y,df_m,optionalParameters,outcome_index){
  run_time<-NA
  YA<- data.frame(Y=Y,A=A)
  TE <- unname(coef(glm(Y ~ A, data = YA))[2])
  
  run_time<-system.time(hima.fit<-hima(Y~A,data.pheno=YA,data.M=df_m,
                                       mediator.type="gaussian",
                                       scale=TRUE,
                                       sigcut=1,
                                       verbose=TRUE,
                                       penalty = "MCP"))[time_index] #penalty default = DBlasso
  
  #HIMA returns NULL - no mediation whatsoever
  if(length(hima.fit$ID)==0) { 
    TIE<-0
    DE<-TE-TIE
    feature_effects<-data.frame("CIE"=rep(0,ncol(df_m)),a=rep(0),b=rep(0),
                                CIE_pvalue=rep(1),
                                CIE_qvalue=rep(1),
                                num_nonzero_feat=rep(0),row.names=colnames(df_m))
  }else{
    
    # TE <- unname(coef(glm(Y ~ A, data = YA))[2])
    TIE<-sum(hima.fit$`alpha*beta`)
    DE<-TE-TIE
    
    feature_effects<-data.frame(row.names=colnames(df_m))
    
    df_hima <- data.frame(
      ID = hima.fit$ID,
      alpha_beta = hima.fit$`alpha*beta`,
      alpha = hima.fit$alpha,
      beta = hima.fit$beta,
      p_value = hima.fit$`p-value`
    )
    rownames(df_hima) <- df_hima$ID
    
    df_hima$CIE_qvalue=p.adjust(df_hima$p_value,method="BH")
    
    feature_effects<-merge(x=feature_effects,y=df_hima,by="row.names",all.x=TRUE)
    row.names(feature_effects)<-feature_effects$Row.names
    feature_effects$Row.names <- NULL
    feature_effects$ID <- NULL
    feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]
    colnames(feature_effects)=c("CIE","a","b","CIE_pvalue","CIE_qvalue")
    feature_effects$num_nonzero_feat <- nrow(df_hima)
    feature_effects$CIE_pvalue[which(is.na(feature_effects$CIE_pvalue))]=1
    feature_effects$CIE_qvalue[which(is.na(feature_effects$CIE_qvalue))]=1
    feature_effects[is.na(feature_effects)] = 0
  }
  return(list(TE=TE,
              DE=DE,
              TIE=TIE,
              run_time=run_time,
              feature_effects=feature_effects,
              TIE_var=NA,
              TIE_pvalue=NA,
              nPC=NA))
}


run_hdma<-function(A,Y,df_m,optionalParameters,outcome_index){
  ## need to standardize microbiome data first
  df_m <- as.data.frame(lapply(df_m, scale))
  
  run_time<-NA
  run_time<-system.time(hdma.fit<-hdma (X=A, Y=Y, M=df_m, COV.XM = NULL, COV.MY = NULL, family = c("gaussian"), method = c("lasso"), topN = NULL,
                                        parallel = FALSE, verbose = TRUE))[time_index] #penalty options are lasso or ridge
  #HDMA returns NULL if no mediation whatsoever
  #"No mediator is identified!"
  if(is.null(hdma.fit)){
    YA<- data.frame(Y=Y,A=A)
    TE <- unname(coef(glm(Y ~ A, data = YA))[2])
    TIE<-0
    DE<-TE-TIE
    feature_effects<-data.frame(row.names=colnames(df_m),"CIE"=rep(0,ncol(df_m)),"a"=rep(0),"b"=rep(0),"CIE_pvalue"=rep(1),
                                CIE_pvalue_Bonferroni=rep(1),CIE_qvalue=rep(1),num_nonzero_feat=rep(0))
  }else{
    #Use output from HDMA function
    TE<-unique(hdma.fit$gamma)
    TIE<-sum(hdma.fit$`alpha*beta`)
    DE<-TE-TIE
    
    hdma.fit$CIE_pvalue_Bonferroni<-p.adjust(hdma.fit$P.value,method="bonferroni",n=nrow(hdma.fit))
    hdma.fit$CIE_qvalue<-p.adjust(hdma.fit$P.value,method="fdr",n=nrow(hdma.fit))
    hdma.fit$num_nonzero_feat<-nrow(hdma.fit)
    
    feature_effects<-data.frame(row.names=colnames(df_m))
    feature_effects<-merge(x=feature_effects,y=hdma.fit[,c("alpha*beta","alpha","beta","P.value","CIE_pvalue_Bonferroni","CIE_qvalue","num_nonzero_feat")],
                           by="row.names",all.x=TRUE)
    row.names(feature_effects)<-feature_effects$Row.names
    feature_effects<-feature_effects[,-1]
    feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]
    colnames(feature_effects)=c("CIE","a","b","CIE_pvalue","CIE_pvalue_Bonferroni","CIE_qvalue","num_nonzero_feat")
    feature_effects$CIE_pvalue[which(is.na(feature_effects$CIE_pvalue))]=1
    feature_effects$CIE_pvalue_Bonferroni[which(is.na(feature_effects$CIE_pvalue_Bonferroni))]=1
    feature_effects$CIE_qvalue[which(is.na(feature_effects$CIE_qvalue))]=1
    feature_effects[is.na(feature_effects)] = 0
  }
  return(list(TE=TE,
              DE=DE,
              TIE=TIE,
              run_time=run_time,
              feature_effects=feature_effects,
              TIE_var=NA,
              TIE_pvalue=NA,
              nPC=NA))
}

#Code for naive method
my_naive<-function(A,Y,df_m,optionalParameters,outcome_index){
  feature_effects<-data.frame("Feature"=colnames(df_m))
  selected_features<-which(apply(df_m,2,sum)!=0)
  
  # scale mediator
  df_m<-as.data.frame(sapply(1:ncol(df_m),function(x) (df_m[,x]-mean(df_m[,x]))/sd(df_m[,x])))
  
  # lm(Mk ~ A), for where Mk is a single microbiome feature
  m1_list<-lapply(selected_features, function(x) lm(df_m[, x] ~ A))
  a_coeff<-do.call(rbind, lapply(names(m1_list), function(x) data.frame(Feature=x,a=m1_list[[x]]$coefficients[2])))
  # lm(Y ~ Mk+A), for where Mk is a single microbiome feature
  m2_list<-lapply(selected_features, function(x) lm(Y~df_m[, x] + A))
  b_coeff<-do.call(rbind, lapply(names(m2_list), function(x) data.frame(Feature=x,b=m2_list[[x]]$coefficients[2])))
  
  feature_effects<-merge(x=feature_effects,a_coeff,by="Feature",all.x=TRUE)
  feature_effects<-merge(x=feature_effects,b_coeff,by="Feature",all.x=TRUE)
  feature_effects[is.na(feature_effects)] = 0
  feature_effects$CIE<-feature_effects$a*feature_effects$b
  row.names(feature_effects)<-feature_effects$Feature
  feature_effects<-feature_effects[,-1]
  feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]
  feature_effects<-feature_effects[,c(3,1,2)]
  TIE<-sum(feature_effects$CIE)
  
  #check on this
  m3<-lm(Y~A)
  TE<-summary(m3)$coefficients["A",1]
  DE<-TE-TIE
  
  permutation_pvalue <-NA
  
  #hypothesis test for naive method 
  if(fPerformHypothesisTesting==TRUE){
    df_a<-NULL
    df_b<-NULL
    set.seed(outcome_index*1000)
    for(i in 1:nIterationsHypothesisTesting){
      #PERFORM SHUFFLING
      
      # #BASE - within features, we shuffle which sample has each value for that feature (shuffle column-wise)
      #Shuffle rows without breaking them
      reordered_rows<-sample(nrow(df_m),replace=FALSE)
      
      # lm(Mk ~ A), for where Mk is a single microbiome feature
      m1_list<-lapply(selected_features, function(x) lm(df_m[reordered_rows, x] ~ A))
      a_coeff_perm<-do.call(rbind, lapply(names(m1_list), function(x) data.frame(Feature=x,a=m1_list[[x]]$coefficients[2])))
      a_coeff_perm$iteration<-i
      df_a<-rbind(df_a,a_coeff_perm)
      
      # lm(Mk ~ A), for where Mk is a single microbiome feature
      m2_list<-lapply(selected_features, function(x) lm(Y~df_m[reordered_rows, x] + A))
      b_coeff_perm<-do.call(rbind, lapply(names(m2_list), function(x) data.frame(Feature=x,b=m2_list[[x]]$coefficients[2])))
      b_coeff_perm$iteration<-i
      df_b<-rbind(df_b,b_coeff_perm)
    }
    
    #combine all data in df_merged
    df_merged<-merge(x=df_a,y=df_b,by=c("Feature","iteration"))  
    rm(a_coeff_perm,b_coeff_perm,df_a,df_b,m1_list,m2_list)
    colnames(df_merged)[3:4]<-c("a_perm","b_perm")
    df_merged<-merge(x=df_merged,y=a_coeff,by=c("Feature"),all.x=TRUE)
    df_merged<-merge(x=df_merged,y=b_coeff,by=c("Feature"),all.x=TRUE)
    
    df_merged$a_perm_b_perm<-df_merged$a_perm*df_merged$b_perm
    df_merged$a_b_perm<-df_merged$a*df_merged$b_perm
    df_merged$a_perm_b<-df_merged$a_perm*df_merged$b
    
    df_test_stats <- merge(
      aggregate(df_merged$a_perm_b_perm, by = list(iteration = df_merged$iteration), FUN = sum),
      aggregate(df_merged$a_perm_b, by = list(iteration = df_merged$iteration), FUN = sum),
      by = "iteration"
    )
    df_test_stats <- merge(
      df_test_stats,
      aggregate(df_merged$a_b_perm, by = list(iteration = df_merged$iteration), FUN = sum),
      by = "iteration"
    )
    colnames(df_test_stats) <- c("iteration", "TIE_a_perm_b_perm", "TIE_a_perm_b", "TIE_a_b_perm")
    
    test_stat <- apply(df_test_stats[,-1], 1, function(row) max(abs(row)))
    
    #The p-value is the number of null test stats that are bigger than the observed value
    #We add a 1 to the numerator and denominator to account for misestimation of the p-value 
    # (for more details see Phipson and Smyth, Permutation P-values should never be zero).
    permutation_pvalue <- (sum(abs(test_stat)>=abs(TIE))+1)/(nIterationsHypothesisTesting+1)
    
    # CIE
    feat_list<-unique(df_merged$Feature)
    
    df_pvalue<-data.frame("Feature"=feat_list)
    df_CIE <- df_merged %>% dplyr::select(Feature,iteration,a_perm_b_perm,a_b_perm,a_perm_b)%>% 
      dplyr::mutate(CIE_test_stats = apply(dplyr::select(., a_perm_b_perm, a_b_perm, a_perm_b), 1, function(x) max(abs(x))))
    df_pvalue$CIE_pvalue <- unlist(lapply(feat_list,function(x) (
      sum(abs(df_CIE$CIE_test_stats[which(df_CIE$Feature==x)]) >= 
            abs(feature_effects$CIE[which(row.names(feature_effects)==x)]))+1)/
        (nIterationsHypothesisTesting+1)))
    
    feature_effects<-merge(x=feature_effects,y=df_pvalue,by.x="row.names",by.y="Feature",all.x=TRUE)
    row.names(feature_effects)<-feature_effects$Row.names
    feature_effects<-feature_effects[,-1]
    feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]
  }else{
    feature_effects$CIE_pvalue<-NA
  }
  return(list(TE=TE,
              DE=DE,
              TIE=TIE,
              run_time=NA,
              feature_effects=feature_effects,
              TIE_var=NA,
              TIE_pvalue=permutation_pvalue,
              nPC=NA))
}


#Wrapper/timing for naive method
run_naive<-function(A,Y,df_m,optionalParameters,outcome_index){
  run_time<-system.time(rslt.naive<-my_naive(A=A,Y=Y,df_m=df_m,optionalParameters=optionalParameters,outcome_index=outcome_index))[time_index]
  rslt.naive$run_time<-run_time
  return(rslt.naive)
}


my_pcr<-function(A,Y,df_m,optionalParameters,outcome_index){
  option<-list(fCenter=TRUE,fScale=TRUE,varThreshold=0.9)
  #if using scaling make sure I remove features that are all 0
  if(option$fScale==TRUE){
    zero_features<-which(apply(df_m,2,sum)==0)
    if(length(zero_features)>0){
      df_m<-df_m[,-zero_features]
    }
  }
  if(fUseNonCompositionalData==TRUE){
    # df_m.trans<-sqrt(df_m)
    df_m.trans<-df_m
  }else{
    df_m.trans<-asin(sqrt(df_m))
  }
  pcs<-prcomp(df_m.trans,center=TRUE,scale=option$fScale)
  cum_prop<-(summary(pcs)$importance)[3,]
  num_pcs<-min(which(cum_prop>option$varThreshold))
  df_z<-as.data.frame(pcs$x[,1:num_pcs])
  colnames(df_z)<-paste0("PC",seq(1,ncol(df_z)))
  
  feat<-1:ncol(df_z)
  names(feat)<-colnames(df_z)
  m1_list<-lapply(feat, function(x) lm(df_z[, x] ~ A))
  a_coeff<-do.call(rbind, lapply(names(m1_list), function(x) data.frame(Feature=x,a=m1_list[[x]]$coefficients[2])))
  m2_list<-lapply(feat, function(x) lm(Y~df_z[, x] + A))
  b_coeff<-do.call(rbind, lapply(names(m2_list), function(x) data.frame(Feature=x,b=m2_list[[x]]$coefficients[2])))
  
  feature_effects<-data.frame("Feature"=colnames(df_z))
  feature_effects<-merge(x=feature_effects,a_coeff,by="Feature",all.x=TRUE)
  feature_effects<-merge(x=feature_effects,b_coeff,by="Feature",all.x=TRUE)
  feature_effects[is.na(feature_effects)] = 0
  feature_effects$CIE<-feature_effects$a*feature_effects$b
  row.names(feature_effects)<-feature_effects$Feature
  feature_effects<-feature_effects[,-1]
  feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]
  feature_effects<-feature_effects[,c(3,1,2)]
  TIE<-sum(feature_effects$CIE)
  m3<-lm(Y~A)
  TE<-summary(m3)$coefficients["A",1]
  DE<-TE-TIE
  
  pcData<-data.frame(row.names=colnames(df_z))
  pcData<-merge(x=pcData,y=t(summary(pcs)$importance),by="row.names",all.x=TRUE)
  pcData<-merge(x=pcData,y=t(pcs$rotation),by.x="Row.names",by.y="row.names",all.x=TRUE)
  row.names(pcData)<-pcData$Row.names
  pcData<-pcData[,-1]
  
  permutation_pvalue<-NA
  
  if(fPerformHypothesisTesting==TRUE){
    # pcr_name<-paste0("scale_",optionalParameters$fScale,"_var_",optionalParameters$varThreshold)
    df_a<-NULL
    df_b<-NULL
    set.seed(outcome_index)
    for(i in 1:nIterationsHypothesisTesting){
      #PERFORM SHUFFLING
      
      # #BASE - within features, we shuffle which sample has each value for that feature (shuffle column-wise)
      #Shuffle rows without breaking them
      df_m_permuted<-df_m[sample(nrow(df_m)),]
      A_permuted<-A
      Y_permuted<-Y
      
      pcs<-prcomp(asin(sqrt(df_m_permuted)),center=TRUE,scale=option$fScale)
      cum_prop<-(summary(pcs)$importance)[3,]
      num_pcs<-min(which(cum_prop>option$varThreshold))
      df_z_permuted<-as.data.frame(pcs$x[,1:num_pcs])
      colnames(df_z_permuted)<-paste0("PC",seq(1,ncol(df_z_permuted)))
      
      # Specify the columns that contain your predictor variables
      feat<-1:ncol(df_z_permuted)
      names(feat)<-colnames(df_z_permuted)
      
      # lm(Mk ~ A), for where Mk is a single microbiome feature
      m1_list<-lapply(feat, function(x) lm(df_z_permuted[, x] ~ A_permuted))
      a_coeff_perm<-do.call(rbind, lapply(names(m1_list), function(x) data.frame(Feature=x,a=m1_list[[x]]$coefficients[2])))
      a_coeff_perm$iteration<-i
      df_a<-rbind(df_a,a_coeff_perm)
      
      # lm(Y ~ Zk+A), for where Zk is a single PC 
      m2_list<-lapply(feat, function(x) lm(Y_permuted~df_z_permuted[, x] + A_permuted))
      b_coeff_perm<-do.call(rbind, lapply(names(m2_list), function(x) data.frame(Feature=x,b=m2_list[[x]]$coefficients[2])))
      b_coeff_perm$iteration<-i
      df_b<-rbind(df_b,b_coeff_perm)
    }
    
    #combine all data in df_merged
    df_merged<-merge(x=df_a,y=df_b,by=c("Feature","iteration"))  
    rm(a_coeff_perm,b_coeff_perm,df_a,df_b,m1_list,m2_list)
    colnames(df_merged)[3:4]<-c("a_perm","b_perm")
    df_merged<-merge(x=df_merged,y=a_coeff,by=c("Feature"),all.x=TRUE)
    df_merged<-merge(x=df_merged,y=b_coeff,by=c("Feature"),all.x=TRUE)
    
    # 3 scenarios
    df_merged$a_perm_b_perm<-df_merged$a_perm*df_merged$b_perm
    df_merged$a_b_perm<-df_merged$a*df_merged$b_perm
    df_merged$a_perm_b<-df_merged$a_perm*df_merged$b
   
    df_test_stats <- merge(
      aggregate(df_merged$a_perm_b_perm, by = list(iteration = df_merged$iteration), FUN = sum),
      aggregate(df_merged$a_perm_b, by = list(iteration = df_merged$iteration), FUN = sum),
      by = "iteration"
    )
    df_test_stats <- merge(
      df_test_stats,
      aggregate(df_merged$a_b_perm, by = list(iteration = df_merged$iteration), FUN = sum),
      by = "iteration"
    )
    colnames(df_test_stats) <- c("iteration", "TIE_a_perm_b_perm", "TIE_a_perm_b", "TIE_a_b_perm")
    
    test_stat <- apply(df_test_stats[,-1], 1, function(row) max(abs(row)))
    
    permutation_pvalue <- (sum(abs(test_stat)>=abs(TIE))+1)/(nIterationsHypothesisTesting+1)
    
  }else{
    # feature_effects$CIE_pvalue<-NA
    permutation_pvalue=NA
  }
  
  return(list(TE=TE,
              DE=DE,
              TIE=TIE,
              run_time=NA,
              feature_effects=feature_effects,
              TIE_var=NA,
              TIE_pvalue=permutation_pvalue,
              nPC=num_pcs,
              pcData=pcData))
}


#As currently working TSS for normalization, AST for transformation
run_pcr<-function(A,Y,df_m,optionalParameters,outcome_index){
  # print(paste0("PCR: fCenter ",optionalParameters$fCenter,", fScale ",optionalParameters$fScale,", varThreshold ",optionalParameters$varThreshold))
  run_time<-system.time(rslt.pcr<-my_pcr(A=A,Y=Y,df_m=df_m,optionalParameters=optionalParameters,outcome_index=outcome_index))[time_index]
  rslt.pcr$run_time<-run_time
  return(rslt.pcr)
}


## parallel computing for method running
runMethod <- function(A,Y,df_m,i,methodName,methodFunction,optionalParameters) {
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  out <- tryCatch(
    {
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      FUN <- match.fun(methodFunction)
      out<-FUN(A=A,
               Y=Y,
               df_m=df_m,
               optionalParameters=optionalParameters,
               outcome_index=i)
      out$runStatus="Success"
      out$method=methodName
      out$outcome_index=i
      out$message=NA
      out
    },
    error=function(e) {
      # message(paste0("Error for index ",i,":"))
      # message(e)
      # Return value in case of error
      return(list(TE=NA,
                  DE=NA,
                  TIE=NA,
                  run_time=NA,
                  feature_effects=NA,
                  TIE_var=NA,
                  TIE_pvalue=NA,
                  nPC=NA,
                  pcData=NA,
                  runStatus="Error",
                  method=methodName,
                  outcome_index=i,
                  message=e))
    }
  )    
  return(out)
}


run_hima_compositional<-function(A,Y,df_m,optionalParameters,outcome_index){
  run_time<-NA
  YA<- data.frame(Y=Y,A=A)
  TE <- unname(coef(glm(Y ~ A, data = YA))[2])
  
  run_time<-system.time(hima.fit<-hima(Y~A,data.pheno=YA,data.M=df_m,
                                       mediator.type="compositional",
                                       sigcut=0.25,
                                       verbose=TRUE))[time_index] #penalty default = DBlasso
  
  #HIMA returns NULL - no mediation whatsoever
  if(length(hima.fit$ID)==0) { 
    TIE<-0
    DE<-TE-TIE
    feature_effects<-data.frame("CIE"=rep(0,ncol(df_m)),a=rep(0),b=rep(0),
                                CIE_pvalue=rep(1),
                                CIE_qvalue=rep(1),
                                num_nonzero_feat=rep(0),row.names=colnames(df_m))
  }else{
    
    # TE <- unname(coef(glm(Y ~ A, data = YA))[2])
    TIE<-sum(hima.fit$`alpha*beta`)
    DE<-TE-TIE
    
    feature_effects<-data.frame(row.names=colnames(df_m))
    
    df_hima <- data.frame(
      ID = hima.fit$ID,
      alpha_beta = hima.fit$`alpha*beta`,
      alpha = hima.fit$alpha,
      beta = hima.fit$beta,
      p_value = hima.fit$`p-value`
    )
    rownames(df_hima) <- df_hima$ID
    df_hima$CIE_qvalue=0.0001 # set a very small q-value for those significant mediators
    # since here use the FDRcut = 0.25
    
    feature_effects<-merge(x=feature_effects,y=df_hima,by="row.names",all.x=TRUE)
    row.names(feature_effects)<-feature_effects$Row.names
    feature_effects$Row.names <- NULL
    feature_effects$ID <- NULL
    feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]
    colnames(feature_effects)=c("CIE","a","b","CIE_pvalue","CIE_qvalue")
    feature_effects$num_nonzero_feat <- nrow(df_hima)
    feature_effects$CIE_pvalue[which(is.na(feature_effects$CIE_pvalue))]=1
    feature_effects$CIE_qvalue[which(is.na(feature_effects$CIE_qvalue))]=1
    feature_effects[is.na(feature_effects)] = 0
  }
  return(list(TE=TE,
              DE=DE,
              TIE=TIE,
              run_time=run_time,
              feature_effects=feature_effects,
              TIE_var=NA,
              TIE_pvalue=NA,
              nPC=NA))
}
# This one automatically switched to DBlasso 



 