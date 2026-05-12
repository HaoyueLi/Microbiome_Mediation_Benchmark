## This script generates child datasets for simulations
# count version for null and spiked 

func_generate_child_data<-function() {
  
  #set working directory
  setwd(dir_wd)
  
  # add libraries
  library("stringr")
  library("parallel")
  
  # dir.create(paste0(dir_output, "data_creation/child_sets/lognormal/"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(dir_output, "data_creation/child_sets/count/"), recursive = TRUE, showWarnings = FALSE)
  # dir.create(paste0(dir_output, "data_creation/child_sets/spiked_features/"), recursive = TRUE, showWarnings = FALSE)
  
  #load simulation data
  load(paste0(dir_output,"data_creation/list_parent.RData"))
  load(paste0(dir_output,"data_creation/sim.RData")) # sim parameters
  
  
  ## Function for varying library sizes given a relative abundance (i.e. NORMALIZED) table and the mean library size
  #function should 1) transpose, 2) normalize, 3) count sampling, 4) transpose, 5) return results to this function here
  funcGetCounts = function(
    # Lognormal draws  table
    mtrxLogNorm,
    seed
  ) {
    # transpose absolute abundance matrix
    # should use p x n matrix
    mtrxLogNorm<-t(mtrxLogNorm) 
    int_number_samples = ncol(mtrxLogNorm)
    MuLibSize = 10.04278
    SDLibSize = 1.112657
    TimesSDIsOutlierLibSize = 3
    # Average read depth
    iReadDepth = 50000
    vdLogMean = MuLibSize
    vdLogSD = SDLibSize
    viThreshold = TimesSDIsOutlierLibSize
    
    #normalize --> relative abundance
    mtrxRelAb = mtrxLogNorm / rep( colSums( mtrxLogNorm ), each = nrow( mtrxLogNorm ) )
    viZero = which( apply( mtrxLogNorm, 2, sum ) == 0 )
    mtrxRelAb[, viZero] = rep( 0, nrow( mtrxRelAb ) )
    
    set.seed(seed)
    # log-normal distribution
    # assume total read counts per sample followed a log-normal distribution
    iLibSize = exp(truncnorm::rtruncnorm(n = int_number_samples,
                                         mean = vdLogMean,
                                         sd = vdLogSD,
                                         b = vdLogMean + viThreshold * vdLogSD)) #b --> upper truncation point
    iLibSize = round(iLibSize / mean(iLibSize) * iReadDepth) #normalized
    mtrxCounts = sapply(1:int_number_samples,
                        function(i) {
                          prob = mtrxRelAb[, i]
                          if(all(prob == 0)) return(rep(0, nrow(mtrxRelAb)))
                          rmultinom(n = 1,
                                    size = iLibSize[i],
                                    prob = prob)
                        })
    colnames(mtrxCounts)<-colnames(mtrxRelAb)
    rownames(mtrxCounts)<-rownames(mtrxRelAb)
    mtrxCounts<-data.frame(t(mtrxCounts))
    return(mtrxCounts)
  }
  
  make_child_data<-function(df_null,df_spiked,feat,param){
    #add null to lognormal matrix
    df_basis<-df_null
    colnames(df_basis)<-paste0("null_",colnames(df_basis))
    #add null to count matrix
    df_count<-funcGetCounts(mtrxLogNorm=df_null,seed=param$dataset_index)
    colnames(df_count)<-paste0("null_",colnames(df_count))
    
    # generate count for spiked version
    for(i in 1:ncol(feat)){
      df_mixed_basis<-cbind(df_null[,-which(colnames(df_null) %in% feat[,i]),drop=F],
                            df_spiked[,which(colnames(df_spiked) %in% feat[,i]),drop=F])
      df_mixed_basis<-df_mixed_basis[,str_order(colnames(df_mixed_basis),numeric=TRUE)]
      colnames(df_mixed_basis)<-paste0(colnames(feat)[i],"_",colnames(df_mixed_basis))
      df_mixed_count<-funcGetCounts(mtrxLogNorm=df_mixed_basis,seed=param$seed)
      
      df_basis<-cbind(df_basis,df_mixed_basis)
      df_count<-cbind(df_count,df_mixed_count)
    }
    return(list(basis=df_basis,count=df_count))
  }
  
  nCores=20
  
  list_child <- mclapply(1:nrow(sim), mc.cores = nCores, function(i) {
    print(paste0("Generate child data: ",i," of ",nrow(sim)))
    
    ## Read parent dataset - get parameters for it
    # index --> which parameter set
    parent_index<-list_parent[[i]]$index
    
    #error checking
    if(i!=parent_index){
      stop("Error in generate_child_data: indexing off")
    }
    
    parent_param<-sim[which(sim$dataset_index==parent_index),]
    
    # get null and spiked parent datasets
    df_parent_null<-list_parent[[i]]$basis_null
    df_parent_spiked<-list_parent[[i]]$basis_spiked
    feat_char<-list_parent[[i]]$feat_char
    
    # identify features to spike
    feat_spiked<-data.frame(spiked=colnames(df_parent_spiked))
    
    # make child bugs
    lChild<-make_child_data(df_null=df_parent_null,
                            df_spiked=df_parent_spiked,
                            feat=feat_spiked,
                            param=parent_param)
    
    return(list(basis=lChild$basis,count=lChild$count,feat_spiked=feat_spiked,feat_char=feat_char))
  })
  
  
  save(list_child, file = paste0(dir_output, "data_creation/list_child.RData"))
  
}
