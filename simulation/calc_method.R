## This script contains all mediation method functions for the project
# no SparseMCMM function (this method was only used for real-world data analysis)

# set working directory
# dir_wd<-paste0(here::here(),"/") #working directory
# setwd(dir_wd)


library("stringr")
library("ccmm")
library("SparseMCMM")
library("readxl")
library("vegan")
source("helper_functions.R")
library("mvtnorm")
library("Matrix")
library("Rfast")
library("matrixStats")
library("ncvreg")
library("HIMA")
library(LDM)


#additional settings
time_index<-1 #which value of system.time to use

### Naive
my_naive<-function(A,Y,df_m,optionalParameters,outcome_index){
  
  feature_effects<-data.frame("Feature"=colnames(df_m))
  
  selected_features<-which(apply(df_m,2,sum)!=0)
  
  my_names<-colnames(df_m)
  df_m<-as.data.frame(sapply(1:ncol(df_m),function(x) (df_m[,x]-mean(df_m[,x]))/sd(df_m[,x])))
  colnames(df_m)<-my_names
  
  ## mediator model
  # lm(Mk ~ A), for where Mk is a single microbiome feature
  m1_list<-lapply(selected_features, function(x) lm(df_m[, x] ~ A))
  a_coeff<-do.call(rbind, lapply(names(m1_list), function(x) data.frame(Feature=x,a=m1_list[[x]]$coefficients[2])))
  
  ## outcome model 
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
  
  # direct effect
  m3<-lm(Y~A)
  TE<-summary(m3)$coefficients["A",1]
  DE<-TE-TIE
  
  permutation_pvalue <-NA
  
  ## hypothesis test for naive method 
  # permutation-based
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
    
    
    ## CIE calculation
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
  }
  
  else{
    ## do not perform hypothesis testing
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


### PCR
my_pcr<-function(A,Y,df_m,optionalParameters,outcome_index){
  
  #if using scaling make sure I remove features that are all 0
  if(optionalParameters$fScale==TRUE){
    zero_features<-which(apply(df_m,2,sum)==0)
    if(length(zero_features)>0){
      df_m<-df_m[,-zero_features]
    }
  }
  if(fUseNonCompositionalData==FALSE & transformationType=="TSS"){
    df_m.trans<-asin(sqrt(df_m))
  }else{
    df_m.trans<-df_m 
  }
  
  pcs<-prcomp(df_m.trans,center=TRUE,scale=optionalParameters$fScale)
  cum_prop<-(summary(pcs)$importance)[3,]
  num_pcs<-min(which(cum_prop>optionalParameters$varThreshold))
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
  
  ## hypothesis testing
  # permutation-based
  if(fPerformHypothesisTesting==TRUE){
    pcr_name<-paste0("scale_",optionalParameters$fScale,"_var_",optionalParameters$varThreshold)
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
      
      pcs<-prcomp(asin(sqrt(df_m_permuted)),center=TRUE,scale=optionalParameters$fScale)
      cum_prop<-(summary(pcs)$importance)[3,]
      num_pcs<-min(which(cum_prop>optionalParameters$varThreshold))
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
    ## no hypothesis testing performed
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
  print(paste0("PCR: fCenter ",optionalParameters$fCenter,", fScale ",optionalParameters$fScale,", varThreshold ",optionalParameters$varThreshold))
  run_time<-system.time(rslt.pcr<-my_pcr(A=A,Y=Y,df_m=df_m,optionalParameters=optionalParameters,outcome_index=outcome_index))[time_index]
  rslt.pcr$run_time<-run_time
  return(rslt.pcr)
}


### HDMA
run_hdma<-function(A,Y,df_m,optionalParameters,outcome_index){
  
  ## need to standardize microbiome data first
  ## for real-world data analysis
  # df_m <- as.data.frame(lapply(df_m, scale))
  
  run_time<-NA
  run_time<-system.time(hdma.fit<-hdma (X=A, Y=Y, M=df_m, COV.XM = NULL, COV.MY = NULL, family = c("gaussian"), method = c("lasso"), topN = NULL,
                                        parallel = FALSE, verbose = TRUE))[time_index] #penalty options are lasso or ridge
  #HDMA returns NULL if no mediation whatsoever
  #"No mediator is identified!"
  if(is.null(hdma.fit)){
    TE<-NA
    TIE<-0
    DE<-NA
    feature_effects<-data.frame(CIE=rep(0,ncol(df_m)),a=rep(0),b=rep(0),
                                CIE_pvalue=rep(1),CIE_pvalue_Bonferroni=rep(1),
                                CIE_qvalue=rep(1),num_nonzero_feat=rep(0),
                                row.names=colnames(df_m))
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


### HIMA Gaussian (non-compositional)
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
    
    #NOTE: not compatible with hypothesis testing
    # hima_features<-as.character(row.names(hima.fit))
    hima_features<-as.character(hima.fit$ID)
    #get different beta values
    YAM<-cbind(Y,A,df_m[,hima_features,drop=F])
    joint_outcome_model<-lm(Y~.,data=YAM)
    DE<-summary(joint_outcome_model)$coefficients["A",1]
    hima.fit$beta<-summary(joint_outcome_model)$coefficients[hima.fit$ID, 1]
    TIE<-sum(hima.fit$`alpha*beta`)
    m3<-lm(Y~A)
    TE<-summary(m3)$coefficients["A",1]
    
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



### CCMM
run_ccmm<-function(A,Y,df_m,optionalParameters,outcome_index){
  M<-as.matrix(df_m)
  
  if(fPerformHypothesisTesting==FALSE){
    run_time<-system.time(rslt.ccmm <- ccmm(y=Y,M=M,tr=A,x=NULL))[time_index]
    feature_effects<-data.frame("CIE"=rslt.ccmm$IDEs,"a"=rep(NA),"b"=rep(NA),
                                "CIE_var"=rep(NA),row.names=colnames(df_m))
    TIE_var<-NA
  }else{
    run_time<-system.time(rslt.ccmm <- ccmm(y=Y,M=M,tr=A,x=NULL,method.est.cov="normal"))[time_index]
    feature_effects<-data.frame("CIE"=rslt.ccmm$IDEs,"a"=rep(NA),"b"=rep(NA),
                                "CIE_var"=(rslt.ccmm$Var.IDEs),row.names=colnames(df_m))
    TIE_var<-unname(rslt.ccmm$Var.TIDE)
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
              TIE_pvalue=NA,
              nPC=NA))
}


### MedTest
run_medtest<-function(A,Y,df_m,optionalParameters,outcome_index){
  set.seed(outcome_index*100)
  
  if(fPerformHypothesisTesting==FALSE){
    TIE_pvalue=NA
    run_time=NA
  }else{
    M<-as.matrix(df_m)
    M.list <- list(BC=as.matrix(vegan::vegdist(M, method="bray")), 
                   JAC=as.matrix(vegan::vegdist(M, 'jaccard', binary=TRUE)),
                   EU=as.matrix(vegan::vegdist(M, method="euclidean")))
    run_time<-system.time(rslt.medtest <- MedOmniTest(A, Y, M.list))[time_index]
    TIE_pvalue<-rslt.medtest$permP
  }
  return(list(TE=NA,
              DE=NA,
              TIE=NA,
              run_time=run_time,
              feature_effects=NA,
              TIE_var=NA,
              TIE_pvalue=TIE_pvalue,
              nPC=NA))
}



### LDM-med
my_ldm <- function(A,Y,df_m, optionalParameters,outcome_index){
  seed <- simOutcome$seed[outcome_index]
  covariates <- data.frame(A=A,Y=Y)
  df_m <- as.data.frame(df_m)
  assign("df_m", df_m, envir = .GlobalEnv)

  ldm.med <- LDM::ldm(formula=df_m~ A+Y, data=covariates,
                      seed=seed, n.cores=1, test.mediation=TRUE,
                      test.omni3=TRUE)

  # extract global p-value
  TIE_pvalue_freq <- ldm.med$med.p.global.freq
  TIE_pvalue_pa <- ldm.med$med.p.global.pa
  TIE_pvalue_tran <- ldm.med$med.p.global.tran
  TIE_pvalue_omni3 <- ldm.med$med.p.global.omni3

  # extract CIE
  feature_effects<-data.frame("Feature"=colnames(df_m),"CIE"=rep(NA),"a"=rep(NA),"b"=rep(NA),row.names=colnames(df_m))
  feature_effects$CIE_qvalue_freq <- ldm.med$med.q.otu.freq
  feature_effects$CIE_qvalue_pa <- ldm.med$med.q.otu.pa
  feature_effects$CIE_qvalue_tran <- ldm.med$med.q.otu.tran
  feature_effects$CIE_qvalue_omni3<- ldm.med$med.q.otu.omni3

  feature_effects <- feature_effects[,-1]
  ## for non-significant features, set q-value to 1
  feature_effects$CIE_qvalue_freq[which(is.na(feature_effects$CIE_qvalue_freq))]=1
  feature_effects$CIE_qvalue_pa[which(is.na(feature_effects$CIE_qvalue_pa))]=1
  feature_effects$CIE_qvalue_tran[which(is.na(feature_effects$CIE_qvalue_tran))]=1
  feature_effects$CIE_qvalue_omni3[which(is.na(feature_effects$CIE_qvalue_omni3))]=1

  feature_effects<-feature_effects[str_order(row.names(feature_effects),decreasing=FALSE,numeric=TRUE),]

  return(list(TE=NA,
              DE=NA,
              TIE=NA,
              run_time=NA,
              feature_effects=feature_effects,
              TIE_var=NA,
              TIE_pvalue_freq=TIE_pvalue_freq,
              TIE_pvalue_pa=TIE_pvalue_pa,
              TIE_pvalue_tran=TIE_pvalue_tran,
              TIE_pvalue_omni3=TIE_pvalue_omni3,
              nPC=NA))
}


run_ldm <- function(A,Y,df_m,optionalParameters,outcome_index) {
  run_time <- system.time(rslt.ldm<-my_ldm(A=A,Y=Y,df_m=df_m,
                                           optionalParameters=optionalParameters,outcome_index=outcome_index))[time_index]

  rslt.ldm$run_time<-run_time
  return(rslt.ldm)
}



### HIMA Compositional
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
  }
  else{
    
    #NOTE: not compatible with hypothesis testing
    hima_features<-as.character(hima.fit$ID)
    #get different beta values
    YAM<-cbind(Y,A,df_m[,hima_features,drop=F])
    joint_outcome_model<-lm(Y~.,data=YAM)
    DE<-summary(joint_outcome_model)$coefficients["A",1]
    hima.fit$beta<-summary(joint_outcome_model)$coefficients[hima.fit$ID, 1]
    TIE<-sum(hima.fit$`alpha*beta`)
    m3<-lm(Y~A)
    TE<-summary(m3)$coefficients["A",1]
    
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



## load parameters
# can also load parameters in this script
# dir_wd <- -paste0(here::here(),"/")
# parameter_file<-"parameters"  #main paramters file
# parms<-read_excel(paste0(dir_wd,parameter_file,".xlsx"),na="NA")
# eq<-paste(paste(parms$parameter_name, parms$parameter_value, sep="<-"), collapse=";")
# eval(parse(text=eq))




