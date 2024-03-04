## COMPARING UKB AND GS ##
#calculate correlation between UKB and GS cohend for global, tier and nodal measures

rm(list=ls())
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr)
library(ggplot2);library(reshape);library(reshape2);library(plyr)
library(ggpubr);library(ggh4x);library(cowplot)

#Import data ====
abb <- readxl::read_excel("nodes_85_numbered.xlsx")

#empirical pval for global, tier and nodal 
load("ukb_faonly_cidi_PT_nodal_pval.RData")
load("stradl_faonly_cidi_GS_PT_nodal_pval.RData")

#hfdr corrected pval for global, tier and nodal 
ukb_pval_hfdr <- readRDS("ukb_nodal_tier_hfdr_pval.rds")
stradl_pval_hfdr <- readRDS("stradl_nodal_tier_hfdr_pval.rds")

#Organise data ====

#compare similarity in results across threshold for each measure 
nm <- c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40")

#rename to make it more interpretable 
org_FUN <- function(nodal_pval, tier_pval, pval_hfdr_ls){
  dat <- lapply(1:7, function(x){
    tmp <- rbind(nodal_pval[[x]], tier_pval[[x]])
    tmp$thresh <- nm[x]
    tmp$measure <- as.character(tmp$dv)
    tmp[grep("^V", tmp$measure),] <- tidyr::extract(tmp[grep("^V", tmp$measure),], col = "measure", into = c(NA,"node_pos",NA,"measure"), "(V)([0-9]{1,2})(\\_[a-z]{3,4}\\_)([a-z]{5}\\_[a-z]{2,3})")
    tmp$measure <- ifelse(tmp$measure == "local_cc", "LCC", 
                                 ifelse(tmp$measure == "local_nef", "LNEF",
                                        ifelse(tmp$measure == "glob_cc", "GCC",
                                                      ifelse(tmp$measure == "glob_eff", "GEFF", ""))))
    
    tmp$measure[grep("^t\\d+_local_nef", tmp$dv)] <- "TIER_LNEF"
    tmp$measure[grep("^t\\d+_local_cc", tmp$dv)] <- "TIER_LCC"
    
    tmp <- merge(abb[2:3], tmp %>% select(-abb), by = "region", all.y=T)
    pval_hfdr <- pval_hfdr_ls[[x]]
    tmp1 <- merge(tmp, pval_hfdr[,c("dv","pval_hfdr","pval_hfdr_sig")], by = "dv") #merge with pval_hfdr 
    return(tmp1)})
  
  names(dat) <- nm
  return(dat)}

ukb_cd <- org_FUN(all_pval_cohend_nodal_cidi_fa, tier_all_pval_cohend_nodal_cidi_fa, ukb_pval_hfdr) 
stradl_cd <- org_FUN(all_pval_cohend_nodal_cidi_GS_fa, tier_all_pval_cohend_nodal_cidi_GS_fa, stradl_pval_hfdr)

#Compare UKB and GS cohen's d for all thresholds (7 thresholds)
compare_cohend_FUN <- function(ukb_dat, stradl_dat){
  
  ukb_cd <- ukb_dat %>% select(abb, measure, thresh, cohen_d)
  colnames(ukb_cd)[4] <- paste0(colnames(ukb_cd)[4], "_ukb")
  stradl_cd <- stradl_dat %>% select(abb, measure, thresh, cohen_d)
  colnames(stradl_cd)[4] <- paste0(colnames(stradl_cd)[4], "_stradl")
  
  #extract common regions
  dat_merge <- merge(ukb_cd, stradl_cd, by = c("abb", "measure", "thresh"))  #%>% 

  plt <- ggplot(dat_merge, aes(x = cohen_d_ukb, y = cohen_d_stradl)) +
    geom_point(aes(color = measure)) +
    ggrepel::geom_text_repel(aes(label = abb, colour = measure), size = 3.5, fontface = "bold") +
    scale_color_brewer("Measure", palette = "Set1") + 
    geom_smooth(method = lm, se = F, color = "black") +
    stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = "top") + 
    labs(x = "UKB Cohen's d", 
         y = "GS Cohen's d", 
         title = paste0("UKB VS GS - Network Density: ", unique(dat_merge$thresh))) + 
    theme(text = element_text(face="bold"),
          axis.text = element_text(face="bold"))
  
  return(plt)
}

ukb_stradl_cohend_fa <- mapply(compare_cohend_FUN, ukb_cd, stradl_cd, SIMPLIFY = F)

