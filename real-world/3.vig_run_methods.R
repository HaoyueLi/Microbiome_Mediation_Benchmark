
## This script runs methods on vignette datasets for the mediation project
# except SparseMCMM and MedTest

# set working directory
# dir_wd<-
# setwd(dir_wd)

library(future)
library(parallel)
library("readxl")
library(dplyr)
library(tidyr)

# load the parameters to run methods
load(paste0(dir_wd, "results/df_vig.RData"))


# Set number of cores
nCores<-30

res_vig <- parallel::mclapply(1:nrow(df_vig), mc.cores = nCores, function(i) {
  source("helper_functions.R")
  source("2.vig_calc_methods_data_hutlab.R")

  message(paste0("Run vignettes: ", i, " of ", nrow(df_vig)))

  params_i <- df_vig[i, ,drop=FALSE]

  if (params_i$method == "ccmm"|params_i$method == "hima_compositional"|params_i$method == "hdma"|params_i$method == "SparseMCMM") {
    my_bugs <- get(load(paste0(dir_wd, params_i$dataset, "/bugs_", params_i$microbiome_feature_level, "_pc.RData")))
  } else {
    my_bugs <- get(load(paste0(dir_wd, params_i$dataset, "/bugs_", params_i$microbiome_feature_level, ".RData")))
  }

  # load(paste0(dir_wd, params_i$dataset, "/metadata.RData"))
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


save(res_vig, file = paste0(dir_wd, "results/res_vig.RData"))




