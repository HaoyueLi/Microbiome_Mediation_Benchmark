## This script generates parent datasets for simulations
# generate binary exposure, null & spiked version of simulated microbiome data
# absolute abundance

func_generate_parent_data<-function(nIterations,
                                    nSample,
                                    nMicrobe,
                                    nPathways,
                                    exposureType,
                                    fZeroInflate,
                                    dirPathways_a,
                                    TIE,
                                    fSimplifyCrossing) {
  
  ## add libraries
  library("sparseDOSSA")
  library("tibble")
  library("tidyr")
  library("magrittr")  
  library("stringr")
  library("dplyr")
  library("readr")
  library("parallel")
  source("helper_functions.R")
 
  nCores=20
  num_tries<-250   #number of times it tries again to create a simulated dataset that will pass QC
  minPercent<-0.1  #Parameter for passing Sparsedossa spike QC check
  
  #set working directory
  setwd(dir_wd)
  
  ## create output folders
  dir.create(dir_output, recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(dir_output, "data_creation"), recursive = TRUE, showWarnings = FALSE)
  
  
  ## create dataset with simulation parameters
  simSetup <- crossing(
    nSample,
    nMicrobe,
    nPathways,
    dirPathways_a,
    fZeroInflate,
    exposureType,
    TIE
  )
  
  # only one parameter varies at a time
  # Figure 4
  # sample size = 100, medium TIE, 10 mediating features, 100% mediated, ++
  if(fSimplifyCrossing==TRUE){
    def_nSample<-100
    def_nPathways<-10
    def_dirPathways_a<-1
    def_TIE<-1
    
    simSetup<-rbind(simSetup[which(simSetup$nSample==def_nSample & simSetup$nPathways==def_nPathways & simSetup$dirPathways_a==def_dirPathways_a & simSetup$TIE==def_TIE),],
                    simSetup[which(simSetup$nSample!=def_nSample & simSetup$nPathways==def_nPathways & simSetup$dirPathways_a==def_dirPathways_a & simSetup$TIE==def_TIE),],
                    simSetup[which(simSetup$nSample==def_nSample & simSetup$nPathways!=def_nPathways & simSetup$dirPathways_a==def_dirPathways_a & simSetup$TIE==def_TIE),],
                    simSetup[which(simSetup$nSample==def_nSample & simSetup$nPathways==def_nPathways & simSetup$dirPathways_a!=def_dirPathways_a & simSetup$TIE==def_TIE),],
                    simSetup[which(simSetup$nSample==def_nSample & simSetup$nPathways==def_nPathways & simSetup$dirPathways_a==def_dirPathways_a & simSetup$TIE!=def_TIE),])
  }
  simSetup$param_set<-seq(1,nrow(simSetup))
  
  ## create iterations for running simulation
  # parameter setup, not real microbiome data
  sim<-as.data.frame(lapply(simSetup,function(x)rep(x,nIterations))) 
  sim<-sim[order(sim$param_set,decreasing=FALSE),] #same param_set stack together
  sim$iteration<-rep(seq(1,nIterations),nrow(simSetup)) #same param_set has n iterations
  rownames(sim)<-seq(1,nrow(sim))
  sim$dataset_index<-seq(1,nrow(sim))
  # save(simSetup, file = paste0(dir_output,"data_creation/simSetup.RData"))
  rm(simSetup)
  
  
  ## generate metadata
  # simulated exposure 
  list_meta<-list()
  sim$standard_deviation_treatment<-0
  for(i in 1:nrow(sim)){
    set.seed(sim$dataset_index[i])
    if(sim$exposureType[i]=="binary"){
      meta<-rbinom(n=sim$nSample[i],size=1,prob=0.5)
    }else if(sim$exposureType[i]=="continuous"){
      # meta<-rnorm(n=sim$nSample[i],mean=1,sd=(1/3))
      meta<-rnorm(n=sim$nSample[i],mean=0,sd=5)
      meta<-(meta-mean(meta))/sd(meta) #normalized
    }else{
      #error
      meta<-rep(NA,sim$nSample[i])
    }
    list_meta[[i]]<-as.matrix(t(meta)) #element of the ith list, 1 by sim$nSample[i] matrix
    sim$standard_deviation_treatment[i]<-ifelse(sd(meta)==0,1,sd(meta))
  }
  
  save(list_meta, file = paste0(dir_output,"data_creation/list_meta.RData"))
  
  sim$seed<-sim$dataset_index*num_tries
  save(sim, file = paste0(dir_output, "data_creation/sim.RData"))
  
  
  ### Generate absolute null & spiked microbiome data
  list_parent <- mclapply(1:nrow(sim), mc.cores = nCores, function(i) {
    source("helper_functions.R")
    print(paste0("Generate parent data: ", i, " of ", nrow(sim)))
    
    fPass <- FALSE
    counter <- 0
    vcErrorMessages <- NULL
    
    while (fPass == FALSE & counter < num_tries) {
      
      #generate null data
      run_null.l <- sparseDOSSA::sparseDOSSA(
        datasetCount = 1,
        number_features = sim$nMicrobe[i],
        number_samples = sim$nSample[i],
        max_percent_outliers = 0,
        UserMetadata = list_meta[[i]],
        percent_spiked = 0,
        spikeStrength = "0",
        spikeCount = "1",
        noZeroInflate = !(sim$fZeroInflate[i]),
        seed = (sim$seed[i] + counter),
        minLevelPercent = 0,
        write_table = FALSE,
        verbose = FALSE
      )
      
      run_null <- run_null.l %>%
        extract_sparseDOSSA(dataset = "basis", set = 1) %>%
        extract2("data") %>% t() %>% data.frame
      
      #check if we have features where we can achieve a sufficient spike strength
      percent_nonzero.vc <- apply(run_null != 0, 2, sum) / nrow(run_null)
      feat_opt <- which(percent_nonzero.vc >= (sim$TIE[i] / sim$nPathways[i]))
      
      #add more criteria here
      null_nonzero <- run_null[, feat_opt, drop = FALSE] != 0 #logical matrix, only contains TRUE / FALSE
      sum_0 <- apply(null_nonzero[which(t(list_meta[[i]]) == 0), , drop = FALSE], 2, sum)
      sum_1 <- apply(null_nonzero[which(t(list_meta[[i]]) == 1), , drop = FALSE], 2, sum)
      sum_0 <- sum_0[which(sum_0 >= sim$nSample[i] * minPercent)]
      sum_1 <- sum_1[which(sum_1 >= sim$nSample[i] * minPercent)]
      feat_opt <- feat_opt[intersect(names(sum_1), names(sum_0))]
      
      if (length(feat_opt) < sim$nPathways[i]) {
        counter <- counter + 1
      } else {
        fPass <- TRUE
        feat <- sort(unname(sample(feat_opt, sim$nPathways[i], replace = FALSE)))
        
        strength_a <- sqrt(sim$TIE[i] / sim$nPathways[i]) / percent_nonzero.vc[feat] * sim$dirPathways_a[i]
        path_vals <- sqrt(sim$TIE[i] / sim$nPathways[i])
        
        #now we calculate delta, which is sparsedossa parameter for spike-in
        #delta=(a*sd(t))/(1-(a*sd(t)))
        strength_delta <- pmin((strength_a * sim$standard_deviation_treatment[i]) /
                                 (1 - strength_a * sim$standard_deviation_treatment[i]), 1000)
        
        # generate both null and spiked version after QC (after feature selection)
        my_spikein.mt <- data.frame(
          feature = feat,
          metadata = rep(1, length(feat)),
          strength = strength_delta
        )
        
        simdata <- sparseDOSSA::sparseDOSSA(
          datasetCount = 1,
          number_features = sim$nMicrobe[i],
          number_samples = sim$nSample[i],
          max_percent_outliers = 0,
          UserMetadata = list_meta[[i]],
          percent_spiked = 0,
          spikeStrength = "0",
          spikeCount = "1",
          noZeroInflate = !(sim$fZeroInflate[i]),
          seed = (sim$seed[i] + counter),
          minLevelPercent = 0,
          write_table = FALSE,
          verbose = FALSE,
          spikein.mt = my_spikein.mt
        )
        
        spiked_features <- paste0("Feature", feat)
        
        basis_null <- simdata %>%
          extract_sparseDOSSA(dataset = "basis", set = 1) %>%
          extract2("data") %>% t() %>% data.frame
        
        basis_spiked <- simdata %>%
          extract_sparseDOSSA(dataset = "basis", set = 3) %>%
          extract2("data") %>% t() %>% data.frame
        
        feat_char <- data.frame(
          feature = spiked_features,
          percent_nonzero = percent_nonzero.vc[feat],
          path_vals = path_vals,
          strength_a = strength_a,
          strength_delta = strength_delta
        )
        
        return_value <- list(
          index = i,
          basis_null = basis_null,
          basis_spiked = basis_spiked[, spiked_features, drop = FALSE],
          feat_char = feat_char
        )
      }
    }
    
    if (fPass == FALSE) {
      return_value <- NULL
      stop("Error in generate_parent_data: simulated data did not pass QC after specified number of iterations")
    }
    
    return(return_value)
  })
  
  save(list_parent, file = paste0(dir_output,"data_creation/list_parent.RData"))
  # each parent is a list
  # each list contains a null and a spiked version
  
}
