
## This script runs Medtest on real-world data
# only for TIE significance

# set working directory
# dir_wd<-
# setwd(dir_wd)

library("stringr")
library("dplyr")
library("ape")
library("readr")
library("SparseMCMM")

outcomes_mlvs<-c("mscore_norm","mscore_norm_residuals")
outcomes_blueberry <- c("cGMP_change_norm")

df_vig<-tidyr::crossing(method=c("medtest"), 
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

###########################
run_medtest<-function(A,Y,df_m,optionalParameters,outcome_index){
  set.seed(outcome_index*100)
  if(fPerformHypothesisTesting==FALSE){
    TIE_pvalue=NA
    run_time=NA
  }else{
    M<-as.matrix(df_m)
    M.list <- list(BC=as.matrix(vegdist(M, method="bray")), 
                   JAC=as.matrix(vegdist(M, 'jaccard', binary=TRUE)),
                   EU=as.matrix(vegdist(M, method="euclidean")))
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

#########################

fPerformHypothesisTesting=TRUE
time_index=1
nCores=3 # can change this number

res_vig <- parallel::mclapply(1:nrow(df_vig), mc.cores = nCores, function(i) {
  source("helper_functions.R")
  
  message(paste0("Run vignettes: ", i, " of ", nrow(df_vig)))
  
  params_i <- df_vig[i, ,drop=FALSE]
  
  my_bugs <- get(load(paste0(dir_wd, params_i$dataset, "/bugs_", params_i$microbiome_feature_level, ".RData")))
  
  meta_obj <- load(paste0(dir_wd, params_i$dataset, "/metadata.RData"))
  df <- get(meta_obj)
  
  rows <- which(!is.na(df[, params_i$outcome]))
  outcome <- df[rows, params_i$outcome]
  treatment <- df[rows, params_i$exposure_var]
  df_m <- my_bugs[rows, ]
  
  run.rslt <- runMethod(
    A = treatment,
    Y = outcome,
    df_m = df_m,
    i = i,
    methodName = params_i$method,
    methodFunction = paste0("run_", params_i$method),
    optionalParameters = NA
  )
  
  
  return(run.rslt)
})

save(res_vig, file = paste0(dir_wd, "results/res_medtest.RData"))

## can directly check the results for TIE p-value


