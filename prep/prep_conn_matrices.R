## PREP MATRICES ##

rm(list=ls())
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr)

#Import matrices in 3D format (85x85xno_of_subj) for cases and controls separately

import_mat_FUN <- function(groupings, mdd = c("cidi","cidi_GS"), 
                           measure = c("fa", "_Raw_FA")){
  
  #form case-control based on mdd grouping 
  grp <- colnames(groupings)[grep(mdd, colnames(groupings))]  
  cols.keep <- c("connectome_id",grp)
  
  cases <- groupings %>% 
    .[.[colnames(.)[grep(grp,colnames(.))]]==1, colnames(.) %in% cols.keep] %>% 
    na.omit() %>% 
    .[order(.$connectome_id, decreasing = F),]
  
  control <- groupings %>% 
    .[.[colnames(.)[grep(grp,colnames(.))]]==0, colnames(.) %in% cols.keep] %>% 
    na.omit() %>% 
    .[order(.$connectome_id, decreasing = F),]

  cases_id <- paste(cases$connectome_id, collapse="|") 
  control_id <- paste(control$connectome_id, collapse="|")
  
  #form group matrices 
  files <-  list.files(getwd(),pattern= paste0("*", measure, ".csv"))
  cases_files <- files[str_detect(files,cases_id)]
  control_files <- files[str_detect(files,control_id)] 
  
  #create 3d matrices for each group - 85nodes*85nodes*subjects
  if(mdd == "cidi"){ #for ukb
    cases_mat <- lapply(cases_files,read.csv,header=F)
    control_mat <- lapply(control_files,read.csv,header=F)
    } 
  else{ #aka for stradl because the files are formatted differently
    cases_mat <- lapply(cases_files,read.csv)
    control_mat <- lapply(control_files,read.csv)
  } 
  
  cases_mat <- array(unlist(cases_mat), dim = c(85,85,length(cases_files)),
                     dimnames = list(c(1:85),c(1:85),c(cases$connectome_id)))
  
  control_mat <- array(unlist(control_mat), dim = c(85,85,length(control_files)),
                       dimnames = list(c(1:85),c(1:85),c(control$connectome_id)))
  
  output <- list(cases_mat,control_mat)
  names(output) <- c("cases_mat","control_mat")
  
  return(output)
  
}

#ukb ====

#mdd groupings
ukb_mdd <- readRDS("ukb_mdd_pheno_compiled.rds")

setwd("C:/Users/gladi/OneDrive/Desktop/connectomes/UKBiobank_connectome_2019-07-29/Raw/")
mat_cidi_fa <- import_mat_FUN(ukb_mdd, mdd = "cidi", measure = "fa") 

save.image("ukb_all_cidi_raw_mat.RData")

#stradl ====

#mdd groupings
stradl_mdd <- readRDS("stradl_mdd_pheno_compiled.rds")
  
setwd("C:/Users/gladi/OneDrive/Desktop/connectomes/STRADL_Structural_Connectome/matrices/")
mat_cidi_GS_fa <- import_mat_FUN(stradl_mdd, mdd = "cidi_GS", measure = "_Raw_FA") 

save.image("stradl_all_cidi_GS_raw_mat.RData")

