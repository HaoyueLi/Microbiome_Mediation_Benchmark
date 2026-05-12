## This script runs simulation functions for the mediation project

#############################
#############################
#USER - Set parameters
dir_wd<-paste0(here::here(),"/") #working directory
parameter_file<-"parameters"  #main paramters file


#############################
#############################
#clear environment
rm(list = setdiff(ls(), c("dir_wd","parameter_file")))

#load libraries
library("readxl")
source("generate_parent_data.R")
source("generate_child_data.R")
source("generate_outcome_data.R")


#set working directory and make folder to store results
setwd(dir_wd)
dir.create(paste0(dir_wd, "sim_data/"), recursive = TRUE, showWarnings = FALSE)

#Extract parameters
parms<-read_excel(paste0(dir_wd,parameter_file,".xlsx"),na="NA")
eq<-paste(paste(parms$parameter_name, parms$parameter_value, sep="<-"), collapse=";")
eval(parse(text=eq))


# RUN SIMULATION
#############################
if(fRunDataGeneration==TRUE){
  func_generate_parent_data(nIterations = nIterations,
                            nSample = nSample,
                            nMicrobe = nMicrobe,
                            nPathways = nPathways,
                            exposureType=exposureType,
                            fZeroInflate=fZeroInflate,
                            dirPathways_a=dirPathways_a,
                            TIE=TIE,
                            fSimplifyCrossing=fSimplifyCrossing)
  func_generate_child_data()
  func_generate_outcome_data(DE=DE,
                             fPerformSampleQC=fPerformSampleQC,
                             dirPathways_b=dirPathways_b,
                             effect_type=effect_type,
                             fUseNonCompositionalData=fUseNonCompositionalData,
                             transformationType=transformationType,
                             interactionsFactor=interactionsFactor,
                             fSimplifyCrossing=fSimplifyCrossing)
}

