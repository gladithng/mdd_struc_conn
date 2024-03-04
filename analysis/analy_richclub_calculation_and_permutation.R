## CALCULATION OF RICH CLUB COEFFICIENTS + PERMUTATION TESTING ##

#measured at group level
#script here for UKB, same for GS (not shown)

rm(list=ls())
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr)
library(brainGraph);library(tnet)

setwd("C:/Users/gladi/OneDrive/Desktop/connectomes/analysis/matrices/")
load("ukb_all_cidi_raw_mat.RData")

#STEP 1: Prepare data
#calculation of rich club coefficient requires group-averaged network for cases and control 

#function to prepare network and apply thresholding
org_FUN <- function(mat, thresh){ 
  #"mat" is a list containing 2 lists - cases_mat and control_mat in 3D format (85x85xno_of_subj)
  
  mat_grp <- apply(mat, c(1,2), mean) 
  diag(mat_grp) <- 0
  
  #apply threshold 
  thresh_FUN <- function(indiv_mat, thresh){  
    
    n <- dim(indiv_mat)[1] #number of nodes
    diag(indiv_mat) <- 0
    
    #since symmetric, ensure symmetry is preserved and half number of removed links 
    indiv_mat[lower.tri(indiv_mat)] <- 0
    ud <- 2
    
    #find all links, sort according to magnitude
    ind <- which(indiv_mat != 0) #get index of all non zero links in matrix
    tb_freq <- as.data.frame(cbind(ind,indiv_mat[ind])) %>% .[order(.$V2, decreasing = T),] 
    num_links <- round((n^2-n)*thresh/ud)
    
    #apply threshold 
    indiv_mat[tb_freq[num_links+1:nrow(tb_freq),1]] <- 0
    indiv_mat <- indiv_mat + t(indiv_mat)
    
    return(indiv_mat)
  }
  mat_mean <- thresh_FUN(mat_grp, 0.35)
  
  #check whether matrices are connected
  ts <- graph_from_adjacency_matrix(mat_mean, weighted = T, mode = "undirected", diag = F)
  check <- components(ts)$no #5
  cat(paste("components: ", check, sep="")) 
  
  #calculate mininum spanning tree and add it back to matrix
  mst <- mst(graph_from_adjacency_matrix(mat_grp, weighted = T, mode = "undirected", diag = F))
  adj_mst <- as.matrix(as_adjacency_matrix(mst,attr = "weight"))
  
  mat_mean[which(adj_mst != 0)] <- adj_mst[which(adj_mst != 0)]
  
  #final check on whether matrices are connected
  ts1 <- graph_from_adjacency_matrix(mat_mean, weighted = T, mode = "undirected", diag = F)
  check <- components(ts1)$no #1
  cat(paste("components: ", check, sep="")) 
  
  return(mat_mean)
}

#threshold matrices using 7 different thresholds  
mat_grp_010 <- lapply(mat_cidi_fa, org_FUN, 0.10)
mat_grp_015 <- lapply(mat_cidi_fa, org_FUN, 0.15)
mat_grp_020 <- lapply(mat_cidi_fa, org_FUN, 0.20)
mat_grp_025 <- lapply(mat_cidi_fa, org_FUN, 0.25)
mat_grp_030 <- lapply(mat_cidi_fa, org_FUN, 0.30)
mat_grp_035 <- lapply(mat_cidi_fa, org_FUN, 0.35)
mat_grp_040 <- lapply(mat_cidi_fa, org_FUN, 0.40)

#STEP 2: Calculate rich club coefs ====
#function here adapted from R tnet package 
#adapted for clarity/own understanding of each step

source("richclub_w_lists_FUN.R") 

#calculate rich club coefficient for cases and controls for each threshold
cidi_fa_rich_ww <- lapply(mget(ls(pattern="mat_grp_*")),
                          richclub_w_lists_FUN,reshuffle="weights",
                          NR=1000,directed=NULL)

setwd("C:/Users/gladi/OneDrive/Desktop/connectomes/analysis/rich_coefs")
save.image("ukb_faonly_cidi_richcoefs.RData")

#STEP 3: Permutation testing ====
#to test the significance of rich club coefficient for each group relative to null 
#to test significance of case-control difference in rich club coefficient at each k 
source("richclub_permutation_lists_FUN.R") 

pval_cidi_fa_ww <- lapply(cidi_fa_rich_ww, richclub_permutation_lists_FUN)

#format output to fit plot function later on
extract_FUN <- function(ls, ls_pos){ #ls_pos refer to number of thresholds you have per group
  tmp <- lapply(1:ls_pos, function(x) lapply(ls,'[',x))
  names(tmp) <- c("0.10","0.15","0.20", "0.25", "0.30", "0.35", "0.40")
  return(tmp)}

ls <- mget(ls(pattern="^pval_cidi_fa")) 
grp_pval_cidi_fa <- extract_FUN(ls,7)
