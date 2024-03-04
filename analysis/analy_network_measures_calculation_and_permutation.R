## CALCULATE GLOBAL, TIER AND NODAL NETWORK MEASURES ## 

#measured at individual level
#script here for UKB, same for GS (not shown)

rm(list=ls()) 
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr);library(plyr)
library(brainGraph);library(tnet); library(igraph);library(pbapply)
library(BiocManager); library(BiocParallel)

#Load data ====

#import raw matrices
load("ukb_all_cidi_raw_mat.RData")

#list of node positions and name 
nodes_85 <- read.csv("fslabels_85node.txt", header=F) %>% rownames_to_column(.)
colnames(nodes_85) <- c("node_pos","region")

#list of nodes and their abbreviations
abb <- readxl::read_excel("nodes_85_numbered.xlsx")

#STEP 1: Calculate global and nodal measures ====
source("nodal_measures_indiv_threshprop_FUN.R") 

ls <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
nodal_cidi_fa <- lapply(1:7, function(x) nodal_measures_indiv_threshprop_FUN(mat_cidi_fa,ls[x], bin = "no")) #approx 1h using bplapply/around 45mins for mclapply on eddie
names(nodal_cidi_fa) <- c("0.10","0.15","0.20", "0.25", "0.30", "0.35","0.40") 

#combine nodal df with covariates and mdd status into one df
#to do for each network measure in each threshold
covs <- readRDS("ukb_covs_conn_only.rds")
cidi_group <- readRDS("ukb_mdd_pheno_compiled.rds") %>% dplyr::select(connectome_id, f.eid, mdd_cidi)

fa_df <- lapply(nodal_cidi_fa, `[[`, 2) %>% #extract list containing all the dataframes
  lapply(., function(dat){
    ls <- Reduce(function(x, y) merge(x, y, by = "ID"), dat)
    ls$group <- ifelse(ls$ID %in% cidi_group$connectome_id[cidi_group$mdd_cidi ==1], "cases", "control")
    ls$group <- as.factor(ls$group)
    return(ls)})

allnodes_cidi_fa <- lapply(fa_df, function(x) merge(x, covs[-2], by.x="ID", by.y="connectome_id"))


#STEP 2: Calculate tier-based measures ====
source("nodal_tier_FUN.R") 

tier_allnodes_cidi_fa <- lapply(allnodes_cidi_fa, function(x) nodal_tier_FUN(x, cohort = "ukb"))

#STEP 3: Permutation testing ====
#applied to each list aka each threshold
#use regression to find standardised beta and then calculate cohen's d 
#derive empirical pvalue here - to apply hFDR later
#outcome: nodal measures; predictors: covariate and mdd status 
library(doParallel);library(parallel);library(doFuture);library(future);library(doRNG)
doFuture::registerDoFuture()
future::plan(multicore, workers = parallel::detectCores()-1)

source("nodal_cohend_wperm_FUN.R") 
ls <- c(0.10,0.15,0.20,0.25,0.30,0.35,0.40) 

#for global and nodal measures
reg_covs <- paste("mri_age","sex",sep = "+")
all_pval_cohend_nodal_cidi_fa <- mapply(nodal_cohend_wperm_FUN, 
                                        allnodes_cidi_fa, ls, 
                                        cohort = "ukb", ls_covs = reg_covs, 
                                        perm_N = 1000, SIMPLIFY = F)

#for tier measures 
tier_all_pval_cohend_nodal_cidi_fa <- mapply(nodal_cohend_wperm_FUN, 
                                             tier_allnodes_cidi_fa, ls, 
                                             cohort = "ukb", ls_covs = reg_covs, 
                                             perm_N = 1000, SIMPLIFY = F)

doParallel::stopImplicitCluster()

#STEP 4: Run hierarchical FDR correction ====
#run hFDR for global, tier and nodal measures 
source("hfdr_FUN.R") 

pval_hfdr <- mapply(hfdr_FUN, allnodes_cidi_fa, all_pval_cohend_nodal_cidi_fa, 
                    tier_all_pval_cohend_nodal_cidi_fa, SIMPLIFY = F)


#STEP 5: Save output for easy loading later ====
save.image("ukb_faonly_cidi_PT_nodal_pval.RData")

