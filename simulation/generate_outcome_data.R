## This script generates outcome data for simulations
# convert into relative abundance (noisy relative abundance after count)
# generates outcome data (continuous)

func_generate_outcome_data<-function(DE,
                                     fPerformSampleQC,
                                     dirPathways_b,
                                     effect_type,
                                     fUseNonCompositionalData,
                                     transformationType,
                                     interactionsFactor,
                                     fSimplifyCrossing) {
  
  ## Additional parameters
  cPseudocount_method<-"feature" #by sample ("sample") or by feature ("feature") 
  # use "feature" since moving things to other category
  
  if(cPseudocount_method!="sample" & cPseudocount_method!="feature"){
    stop("Error: Pseudocount must be performed either by sample or by feature")
  }
  
  ## set working directory
  setwd(dir_wd)
  
  ## add libraries
  library("tidyr")
  library("parallel")
  library("boot")
  library("dplyr")
  source("helper_functions.R")
  
  nCores <- 20
  
  dir.create(paste0(dir_output, "main/"), recursive = TRUE, showWarnings = FALSE)
  
  ## load simulation data
  load(paste0(dir_output, "data_creation/list_child.RData")) # count version
  load(paste0(dir_output, "data_creation/list_meta.RData")) # exposure
  load(paste0(dir_output, "data_creation/sim.RData")) # simulation parameter sets
  
  #add additional parameters for making the output
  spiked_feature_type<-colnames(list_child[[1]]$feat_spiked)
  
  simSetupOutcome<-crossing(effect_type,
                            spiked_feature_type,
                            DE,
                            dirPathways_b)
  
  
  if(fSimplifyCrossing==TRUE){
    simSetupOutcome<-simSetupOutcome[which(simSetupOutcome$DE==0|simSetupOutcome$dirPathways_b==1),]
    def_nSample<-100
    def_nPathways<-10
    def_dirPathways_a<-1
    def_TIE<-1
    
    # base scenario
    my_rows<-which(sim$nSample==def_nSample & sim$nPathways==def_nPathways & sim$dirPathways_a==def_dirPathways_a & sim$TIE==def_TIE)
    simOutcome<-cbind(sim[my_rows[1],],simSetupOutcome,row.names = NULL)
    # all iterations
    if(length(my_rows)>1){
      for(i in 2:length(my_rows)){
        simOutcome<-rbind(simOutcome,cbind(sim[my_rows[i],],simSetupOutcome,row.names = NULL))
      }
    }
    # these include varying DE and direction of b 
    
    # varying one parameter case
    sim<-sim[-my_rows,]
    if(nrow(sim)>0){
      
      # vary direction of a 
      my_rows<-which(sim$nSample==def_nSample & sim$nPathways==def_nPathways & sim$dirPathways_a==-1 & sim$TIE==def_TIE)
      simOutcome<-rbind(simOutcome,cbind(sim[my_rows,],effect_type="spiked","spiked_feature_type"="spiked",DE=0,dirPathways_b=-1))
      
      # vary sample size & number of mediators & TIE
      sim$effect_type<-"spiked"
      sim$spiked_feature_type<-"spiked"
      sim$DE<-0
      sim$dirPathways_b<-1
      simOutcome<-rbind(simOutcome,sim)
    }
    
    simOutcome$param_set_new<-case_when(simOutcome$DE==1 ~2,
                                        simOutcome$DE==3 ~3,
                                        simOutcome$dirPathways_a==-1 & simOutcome$dirPathways_b==1 ~4,
                                        simOutcome$dirPathways_a==1 & simOutcome$dirPathways_b==-1 ~5,
                                        simOutcome$dirPathways_a==-1 & simOutcome$dirPathways_b==-1~6,
                                        simOutcome$nSample==50~7,
                                        simOutcome$nSample==200 ~8,
                                        simOutcome$nSample==400 ~9,
                                        simOutcome$nPathways==5 ~10,
                                        simOutcome$nPathways==20 ~11,
                                        simOutcome$TIE==0.5 ~12,
                                        simOutcome$TIE==1.5~13,
                                        TRUE~1)
  }else{
    #add additional parameters for making the output
    simOutcome<-cbind(sim[1,],simSetupOutcome,row.names = NULL)
    if(nrow(sim)>1){
      for(i in 2:nrow(sim)){
        simOutcome<-rbind(simOutcome,cbind(sim[i,],simSetupOutcome,row.names = NULL))
      }
    }
    simOutcome$param_set_new<-NA
  }
  
  ## calculate derived variables
  simOutcome$TIE<-simOutcome$TIE*simOutcome$dirPathways_a*simOutcome$dirPathways_b
  simOutcome$TIE<-ifelse(simOutcome$effect_type=="spiked",simOutcome$TIE,0)
  simOutcome$TE<-simOutcome$TIE+simOutcome$DE
  simOutcome$outcome_index<-seq(1,nrow(simOutcome))
  simOutcome$dirTIE<-simOutcome$dirPathways_a*simOutcome$dirPathways_b
  simOutcome$empTIE_bin<-0
  simOutcome$empTIE_lm<-0
  
  
  ## make list of spiked features that will be used for looking at CIEs
  simOutcome_features<-parallel::mclapply(1:nrow(simOutcome),mc.cores = nCores,function(i) {
    print(paste0("Generate outcome data: feature list ",i," of ",nrow(simOutcome)))
    dataset_index<-simOutcome$dataset_index[i]
    feat<-list_child[[dataset_index]][["feat_spiked"]]
    
    return(list(feat=as.character(feat[,simOutcome$spiked_feature_type[i]]),feat_char=list_child[[dataset_index]]$feat_char))
  })
  
  save(simOutcome_features, file = paste0(dir_output, "main/simOutcome_features.RData"))
  
  
  func_add_pseudocount<-function(df){
    if(cPseudocount_method=="sample"){
      #by sample pseudocount
      df[which(apply(df,1,sum)==0),]<-1
      for(i in 1:nrow(df)){
        row_i<-df[i,]
        min_val<-min(row_i[row_i>0])
        zero_vals<-which(row_i==0)
        if(length(zero_vals)>0){
          df[i,zero_vals]<-min_val/2
        }
      }
    }else if(cPseudocount_method=="feature"){
      #by feature pseudocount
      df[,which(apply(df,2,sum)==0)]<-1
      for(i in 1:ncol(df)){
        col_i<-df[,i]
        min_val<-min(col_i[col_i>0])
        zero_vals<-which(col_i==0)
        if(length(zero_vals)>0){
          df[zero_vals,i]<-min_val/2
        }
      }
    }else{
      stop("Error: Pseudocount must be performed either by sample or by feature")
    }
    return(df)
  }
  
  
  # Abundance filtering - filter out bugs (columns) that are not >0.1% in 5% of samples
  abund_filter<-function(bugs_count){
    bugs<-bugs_count/apply(bugs_count,1,sum)
    col_keep_sum<-rep(NA,ncol(bugs))
    col_keep<-rep(TRUE,ncol(bugs))
    for(i in 1:ncol(bugs)){
      bugsum<-sum(bugs[,i]>=(0.1/100),na.rm=TRUE)
      #must be in at least 5% of samples
      if(bugsum<floor(nrow(bugs)*0.05)){
        col_keep[i]<-FALSE
      }
    }
    other_bugs<-!col_keep
    other<-rep(0,nrow(bugs))
    
    for(i in 1:ncol(bugs_count)){
      if(other_bugs[i]){
        other<-other+bugs_count[,i]
      }
    }
    bugs_count<-bugs_count[,col_keep]
    bugs_count<-cbind(bugs_count,other)
    
    return(bugs_count)
  }
  
  
  ## make relative abundance data that (some) methods will be run on
  func_make_bugs_relab<-function(index,fPseudocount){
    dataset_index<-simOutcome$dataset_index[index]
    df_count<-list_child[[dataset_index]]$count
    #test whether to use null matrix or not
    if(simOutcome$effect_type[index]=="null_ab"|simOutcome$effect_type[index]=="null_a"){
      df<-df_count[,grep("^null_",colnames(df_count))]
      #if null_b, or spiked then use the count features that correspond with the spiked_feature_type
    } else{
      df<-df_count[,grep(paste0("^",simOutcome$spiked_feature_type[index],"_"),colnames(df_count))]
    }
    colnames(df)<-vapply(strsplit(colnames(df), "_", fixed = TRUE), "[", "", 2)
    
    if(fPerformSampleQC==TRUE){
      df<-abund_filter(bugs_count=df)
    }
    
    #check whether to add in a pseudocount and add in pseudocount if necessary > shouldn't interact with  
    #prevalence filtering since we're applying pseudocount per feature
    if(fPseudocount==TRUE){
      df<-func_add_pseudocount(df=df)
      # df[df == 0] <- 1
    }
    
    # different transformations to relative abundance
    # only use TSS in this project
    if(transformationType=="TSS"){
      df_relab<-TSSnorm(features=df)
    }else if(transformationType=="CSS"){
      df_relab<-CSSnorm(features=df)
    }else if(transformationType=="CLR"){
      df_relab<-CLRnorm(features=(TSSnorm(features=df)))
    }else{
      df_relab<-TSSnorm(features=df)
    }
    return(df_relab)
  }
  
  
  ## Generate abundance without PC
  if(fUseNonCompositionalData==TRUE){
    
    ## make mediator data that is not compositional (based on the basis)
    # absolute abundance with no pseudocount
    simOutcome_bugs_relab<-parallel::mclapply(1:nrow(simOutcome),mc.cores=nCores,function(i) {
      
      print(paste0("Generate outcome data: bugs non-compositional ",i," of ",nrow(simOutcome)))
      dataset_index<-simOutcome$dataset_index[i]
      df_basis<-list_child[[dataset_index]]$basis
      feat<-list_child[[dataset_index]]$feat_spiked
      
      basis_null<-df_basis[,grep("^null_",colnames(df_basis))]
      basis_spiked<-df_basis[,grep(paste0("^",simOutcome$spiked_feature_type[i],"_"),colnames(df_basis))]
      sdM0<-vapply(basis_null, function(x) sd(x[x!=0]), numeric(1))
      sdM0<-ifelse(sdM0==0,1,sdM0)
      sdM0<-ifelse(is.na(sdM0),1,sdM0)
      names(sdM0)<-vapply(strsplit(names(sdM0), "_", fixed = TRUE), "[", "", 2)
      
      #if null_a or null_ab then use null scaled basis matrix
      if(simOutcome$effect_type[i]=="null_a"|simOutcome$effect_type[i]=="null_ab"){
        basis_null_scaled<-data.frame(t(t(basis_null)/sdM0))
        colnames(basis_null_scaled)<-vapply(strsplit(colnames(basis_null_scaled), "_", fixed = TRUE), "[", "", 2)
        
        if (fPerformSampleQC == TRUE) {
          basis_null_scaled <- abund_filter(basis_null_scaled)
        }
        return(basis_null_scaled)
      }
      
      #if null_b, or spiked then use the spiked scaled basis matrix
      else{
        basis_spiked_scaled<-data.frame(t(t(basis_spiked)/sdM0))
        colnames(basis_spiked_scaled)<-vapply(strsplit(colnames(basis_spiked_scaled), "_", fixed = TRUE), "[", "", 2)
        if (fPerformSampleQC == TRUE) {
          basis_spiked_scaled <- abund_filter(basis_spiked_scaled)
        }
        return(basis_spiked_scaled)
      }
    })
    
    save(simOutcome_bugs_relab, file = paste0(dir_output, "main/simOutcome_bugs_noncomp_noPC.RData"))
    
    #this error arises due to the parallel workers timing out
    if(length(which(lengths(simOutcome_bugs_relab)==0))!=0){
      stop("Error in generate outcome data: ",length(which(lengths(simOutcome_bugs_relab)==0))," relative abundances are null.")
    } 
    
    rm(simOutcome_bugs_relab)
  }
  else{
    # Generate relative abundance without pseudocount
    # compositional data
    simOutcome_bugs_relab<-parallel::mclapply(1:nrow(simOutcome),mc.cores=nCores,function(i) {
      print(paste0("Generate outcome data: bug relab (no PC) ",i," of ",nrow(simOutcome)))
      retval<-func_make_bugs_relab(index=i,fPseudocount=FALSE)
      return(retval)
    })

    save(simOutcome_bugs_relab, file = paste0(dir_output, "main/simOutcome_bugs_relab_noPC.RData"))
    
    # this error arises due to the parallel workers timing out
    if(length(which(lengths(simOutcome_bugs_relab)==0))!=0){
      stop("Error in generate outcome data: ",length(which(lengths(simOutcome_bugs_relab)==0))," relative abundances are null.")
    }
    rm(simOutcome_bugs_relab)
  
  }
  
  
  ## Generate abundance with PC
  if(fUseNonCompositionalData==TRUE) {
    
    ## for non-compositional data
    # based on basis (absolute abundance)
    simOutcome_bugs_relab<-parallel::mclapply(1:nrow(simOutcome),mc.cores=nCores,function(i) {
      
      print(paste0("Generate outcome data: bugs non-compositional (with PC) ",i," of ",nrow(simOutcome)))
      dataset_index<-simOutcome$dataset_index[i]
      df_basis<-list_child[[dataset_index]]$basis
      feat<-list_child[[dataset_index]]$feat_spiked
      
      basis_null<-df_basis[,grep("^null_",colnames(df_basis))]
      basis_spiked<-df_basis[,grep(paste0("^",simOutcome$spiked_feature_type[i],"_"),colnames(df_basis))]
      sdM0<-vapply(basis_null, function(x) sd(x[x!=0]), numeric(1))
      sdM0<-ifelse(sdM0==0,1,sdM0)
      sdM0<-ifelse(is.na(sdM0),1,sdM0)
      names(sdM0)<-vapply(strsplit(names(sdM0), "_", fixed = TRUE), "[", "", 2)
      
      #if null_a or null_ab then use null scaled basis matrix
      if(simOutcome$effect_type[i]=="null_a"|simOutcome$effect_type[i]=="null_ab"){
        basis_null_scaled<-data.frame(t(t(basis_null)/sdM0))
        colnames(basis_null_scaled)<-vapply(strsplit(colnames(basis_null_scaled), "_", fixed = TRUE), "[", "", 2)
        if (fPerformSampleQC == TRUE) {
          basis_null_scaled <- abund_filter(basis_null_scaled)
        }
      
        basis_null_scaled<-func_add_pseudocount(df=basis_null_scaled)

        return(basis_null_scaled)
      }
      
      #if null_b, or spiked then use the spiked scaled basis matrix
      else{
        basis_spiked_scaled<-data.frame(t(t(basis_spiked)/sdM0))
        colnames(basis_spiked_scaled)<-vapply(strsplit(colnames(basis_spiked_scaled), "_", fixed = TRUE), "[", "", 2)
        if (fPerformSampleQC == TRUE) {
          basis_spiked_scaled <- abund_filter(basis_spiked_scaled)
        }
        
        basis_spiked_scaled<-func_add_pseudocount(df=basis_spiked_scaled)
        
        return(basis_spiked_scaled)
      }
    })
    
    save(simOutcome_bugs_relab, file = paste0(dir_output, "main/simOutcome_bugs_noncomp_PC.RData"))
    
    #this error arises due to the parallel workers timing out
    if(length(which(lengths(simOutcome_bugs_relab)==0))!=0){
      stop("Error in generate outcome data: ",length(which(lengths(simOutcome_bugs_relab)==0))," relative abundances are null.")
    } 
    rm(simOutcome_bugs_relab)
  }
  else {
    # Generate relative abundance with pseudocount
    # compositional data
    simOutcome_bugs_relab<-parallel::mclapply(1:nrow(simOutcome),mc.cores=nCores,function(i) {
      print(paste0("Generate outcome data: bug relab (with PC) ",i," of ",nrow(simOutcome)))
      retval<-func_make_bugs_relab(index=i,fPseudocount=TRUE)
      return(retval)
    })
    
    save(simOutcome_bugs_relab, file = paste0(dir_output, "main/simOutcome_bugs_relab_PC.RData"))
    
    #this error arises due to the parallel workers timing out
    if(length(which(lengths(simOutcome_bugs_relab)==0))!=0){
      stop("Error in generate outcome data: ",length(which(lengths(simOutcome_bugs_relab)==0))," relative abundances are null.")
    }
    rm(simOutcome_bugs_relab)
    
  }
  
  
  ## Get metrics on the spiked features (run on non-PC corrected since that is the "true" data)
  if(fUseNonCompositionalData==FALSE){
    
    load(file = paste0(dir_output, "main/simOutcome_bugs_relab_noPC.RData"))
    simOutcome_feature_metrics<-parallel::mclapply(1:nrow(simOutcome),mc.cores=nCores,function(i) {
      print(paste0("Generate outcome data: feature metrics ",i," of ",nrow(simOutcome)))
      df_relab<-simOutcome_bugs_relab[[i]]
      #subset to spiked features that remain after filtering
      sfeat_remaining<-simOutcome_features[[i]]$feat[which(simOutcome_features[[i]]$feat %in% colnames(df_relab))]
      if(length(sfeat_remaining)>0){
        df_relab<-df_relab[,sfeat_remaining,drop=F]
        ret_value<-data.frame("mean"=apply(df_relab,2,mean),
                              "sd"=apply(df_relab,2,sd),
                              "percent_nonzero"=apply(df_relab!=0,2,sum)/nrow(df_relab))
      }else{
        ret_value<-data.frame("mean"=NA,
                              "sd"=NA,
                              "percent_nonzero"=NA)
      }
      return(ret_value)
    })
    
    save(simOutcome_feature_metrics, file = paste0(dir_output, "main/simOutcome_feature_metrics.RData"))
  }
  rm(simOutcome_bugs_relab,simOutcome_features,simOutcome_feature_metrics,inherits=TRUE)
  
  
  ## make outcome (Y) data, and save off exposure (A) data
  simOutcome_covariates<-parallel::mclapply(1:nrow(simOutcome),mc.cores=nCores,function(i) {
    print(paste0("Generate outcome data: covariates ",i," of ",nrow(simOutcome)))
    
    dataset_index<-simOutcome$dataset_index[i]
    df_basis<-list_child[[dataset_index]]$basis
    feat<-list_child[[dataset_index]]$feat_spiked
    strength_b<-(list_child[[dataset_index]]$feat_char$path_vals)*simOutcome$dirPathways_b[i]
    
    set.seed(simOutcome$dataset_index[i])
    error<-rnorm(simOutcome$nSample[i],mean=0,sd=1)
    
    exposure<-t(list_meta[[dataset_index]])
    
    # if null_ab or null_b, only direct effect from treatment
    if(simOutcome$effect_type[i]=="null_ab"|simOutcome$effect_type[i]=="null_b"){
      outcome=((exposure*simOutcome$DE[i])+error)
    }
    
    # if null_a, or spiked then use the absolute abundance that correspond with the effect_type
    else{
      basis_null<-df_basis[,grep("^null_",colnames(df_basis))]
      basis_spiked<-df_basis[,grep(paste0("^",simOutcome$spiked_feature_type[i],"_"),colnames(df_basis))]
      sdM0<-vapply(basis_null, function(x) sd(x[x!=0]), numeric(1))
      sdM0<-ifelse(sdM0==0,1,sdM0)
      sdM0<-ifelse(is.na(sdM0),1,sdM0)
      names(sdM0)<-vapply(strsplit(names(sdM0), "_", fixed = TRUE), "[", "", 2)
      feat_nonzero<-feat[,simOutcome$spiked_feature_type[i]]
      if(sum(is.na(sdM0[which(names(sdM0) %in% feat_nonzero)]))>0){
        stop("Error in generate outcome data: cannot calculate standard deviation of spiked feature")
      }
      if(simOutcome$effect_type[i]=="spiked"){
        basis_spiked_scaled<-data.frame(t(t(basis_spiked)/sdM0))
        colnames(basis_spiked_scaled)<-vapply(strsplit(colnames(basis_spiked_scaled), "_", fixed = TRUE), "[", "", 2)
        
        sum_interaction<-rep(0,nrow(basis_spiked_scaled))
        ## add interactions between microbiome
        # currently not used
        if(interactionsFactor>0){
          n_interactions<-(length(feat_nonzero)/5)*interactionsFactor
          pairs<-combn(feat_nonzero, 2)
          pairs<-pairs[,sample(seq(1:ncol(pairs)),size=n_interactions,replace=FALSE),drop=F]
          names(strength_b)<-feat_nonzero
          for(j in 1:ncol(pairs)){
            
            sub<-basis_spiked_scaled[,which(colnames(basis_spiked_scaled) %in% pairs[,j]),drop=FALSE]
            
            sum_interaction<-sum_interaction+mean(strength_b)*sub[,1]*sub[,2]
          }
        }
        outcome=((exposure*simOutcome$DE[i])+error+
                   apply(as.data.frame(t(t(basis_spiked_scaled[,which(colnames(basis_spiked_scaled) %in% feat_nonzero),drop=FALSE])*strength_b)),1,sum)+
                   sum_interaction)
      }else{
        #simOutcome$effect_type[i]=="null_a" --> use null features
        # since no spiked in association between treatment and features
        basis_null_scaled<-data.frame(t(t(basis_null)/sdM0))
        colnames(basis_null_scaled)<-vapply(strsplit(colnames(basis_null_scaled), "_", fixed = TRUE), "[", "", 2)
        outcome=((exposure*simOutcome$DE[i])+error+
                   apply(as.data.frame(t(t(basis_null_scaled[,which(colnames(basis_null_scaled) %in% feat_nonzero),drop=FALSE])*strength_b)),1,sum))
      }
    }
    df<-data.frame(exposure=exposure,outcome=outcome)
    
    return(df)
    
  })
  
  
  ## make outcome (Y) data, and save off exposure (A) data
  for(i in 1:nrow(simOutcome)){
    df<-simOutcome_covariates[[i]]
    simOutcome$empTIE_bin[i]<-mean(df$outcome[which(df$exposure==1)])-mean(df$outcome[which(df$exposure==0)])-simOutcome$DE[i]
    simOutcome$empTIE_lm[i]<-coef(lm(df$outcome~df$exposure))[2]-simOutcome$DE[i]
  }
  save(simOutcome_covariates, file = paste0(dir_output, "main/simOutcome_covariates.RData"))
  save(simOutcome, file = paste0(dir_output, "data_creation/simOutcome.RData"))
  

}
