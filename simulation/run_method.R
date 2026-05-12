
## This script contains the function to run mediation functions for simulations
# use parallel computing

dir_wd <- paste0(here::here(),"/") #working directory
setwd(dir_wd)
dir_output <- paste0(dir_wd,"/sim_data/")

library("dplyr")
library("parallel")


runMethod <- function(i,methodName,methodFunction,optionalParameters) {
  
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  out <- tryCatch(
    {
      # 'tryCatch()' will return the last evaluated expression 
      # in case the "try" part was completed successfully
      FUN <- match.fun(methodFunction)
      out<-FUN(A=simOutcome_covariates[[i]]$exposure,
               Y=simOutcome_covariates[[i]]$outcome,
               df_m=simOutcome_bugs_relab[[i]],
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


runParallel <- function(methodName, methodFunction, optionalParameters) {
  library(parallel)
  
  method_output_dir <- file.path(dir_output, "results_metrics", methodName)
  if (!dir.exists(method_output_dir)) {
    dir.create(method_output_dir, recursive = TRUE)
  }
  
  load(paste0(dir_output, "main/simOutcome_covariates.RData")) # simOutcome_covariates
  load(paste0(dir_output, "data_creation/simOutcome.RData"))   # simOutcome
  
  # can adjust the number of cores
  results_metrics <- parallel::mclapply(1:nrow(simOutcome), mc.cores = 40, function(i) { 
    data.table::setDTthreads(1)
    
    result_file <- file.path(method_output_dir, paste0("result_", i, ".RData"))
    if (file.exists(result_file)) {
      message(paste("Skipping outcome", i, "— result already exists."))
      return(NULL)
    }
    
    source("helper_functions.R")
    source("calc_method.R")
    
    if ((simOutcome$exposureType[i] == "continuous" & dfMethods[methodName, "fContinuousExposure"] == FALSE) |
        (simOutcome$exposureType[i] == "binary" & dfMethods[methodName, "fBinaryExposure"] == FALSE)) {
      res <- list(
        TE = NA, DE = NA, TIE = NA, run_time = NA, feature_effects = NA,
        TIE_var = NA, TIE_pvalue = NA, nPC = NA, pcData = NA,
        runStatus = "Error", method = methodName, outcome_index = i,
        message = paste0(methodName, " can't be run for exposure type ", simOutcome$exposureType[i])
      )
    } else if (simOutcome$nMicrobe[i] > dfMethods[methodName, "nMaxFeatures"]) {
      res <- list(
        TE = NA, DE = NA, TIE = NA, run_time = NA, feature_effects = NA,
        TIE_var = NA, TIE_pvalue = NA, nPC = NA, pcData = NA,
        runStatus = "Error", method = methodName, outcome_index = i,
        message = paste0(methodName, " isn't run for this many features: ", simOutcome$nMicrobe[i])
      )
    } else {
      res <- runMethod(i = i, methodName = methodName, 
                       methodFunction = methodFunction, 
                       optionalParameters = optionalParameters)
    }
    
    save(res, file = file.path(method_output_dir, paste0("result_", i, ".RData")))
    rm(res)
    gc()
    return(NULL)
  })
  
  return(results_metrics)
}


## load methods
dfMethods<-as.data.frame(readxl::read_excel("methods.xlsx",na="NA"))
row.names(dfMethods)<-dfMethods$name
## can further select the methods to use
# dfMethods <- dfMethods %>% filter(name %in% c("pcr","hdma","naive"))

## load metadata
load(paste0(dir_output, "main/simOutcome_covariates.RData")) # simOutcome_covariates
load(paste0(dir_output, "data_creation/simOutcome.RData")) # simOutcome

## load parameters
parameter_file<-"parameters"  #main paramters file
parms<-read_excel(paste0(dir_wd,parameter_file,".xlsx"),na="NA")
eq<-paste(paste(parms$parameter_name, parms$parameter_value, sep="<-"), collapse=";")
eval(parse(text=eq))


## parallel computing
for(j in 1:nrow(dfMethods)){
  start.time <- Sys.time()
  rm(simOutcome_bugs_relab)
  
  if(dfMethods$fPseudocount[j]==TRUE){
    load(paste0(dir_output, "main/simOutcome_bugs_relab_PC.RData"))
  }else{
    load(paste0(dir_output, "main/simOutcome_bugs_relab_noPC.RData"))
  }

  runParallel(methodName=dfMethods$name[j],
              methodFunction=dfMethods$functionName[j],
              optionalParameters=eval(parse(text=paste('list(', dfMethods$optionalParameters[j], ')'))))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0("Time to run ",dfMethods$name[j],": ",round(time.taken,3)))
}


