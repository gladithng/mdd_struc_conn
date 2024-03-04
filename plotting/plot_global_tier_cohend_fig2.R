## PLOT BETA VALUES FOR GLOBAL AND TIER MEASURES ## 
#only for best performing threshold aka 35%

rm(list=ls())
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr)
library(ggplot2);library(reshape);library(reshape2);library(plyr)
library(ggpubr);library(ggh4x);library(cowplot)

#Import data ====
abb <- readxl::read_excel("nodes_85_numbered.xlsx")

#empirical pval
load("ukb_faonly_cidi_PT_nodal_pval.RData")
load("stradl_faonly_cidi_GS_PT_nodal_pval.RData")

#hfdr corrected pval
ukb_pval_hfdr <- readRDS("ukb_nodal_tier_hfdr_pval.rds")
stradl_pval_hfdr <- readRDS("stradl_nodal_tier_hfdr_pval.rds")

#node labels, abbreviations
nodes_85 <- read.csv("fslabels_85node.txt", header=F) %>% rownames_to_column(.)
colnames(nodes_85) <- c("node_pos","region")

abb <- readxl::read_excel("nodes_85_numbered.xlsx")

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


#Plot global cohend  ====
best_thresh_FUN <- function(ukb_dat, stradl_dat){
  
  #reorg dat for plotting purposes 
  ukb_dat$abb[grep("glob_cc", ukb_dat$dv)] <- "GCC"
  ukb_dat$abb[grep("glob_eff", ukb_dat$dv)] <- "GEFF"
  ukb_dat$measure <- as.character(ukb_dat$measure)
  ukb_dat$measure[grep("glob_cc|glob_eff", ukb_dat$dv)] <- "GLOBAL"
  ukb_dat$cohort <- "UKB"
  
  stradl_dat$abb[grep("glob_cc", stradl_dat$dv)] <- "GCC"
  stradl_dat$abb[grep("glob_eff", stradl_dat$dv)] <- "GEFF"
  stradl_dat$measure <- as.character(stradl_dat$measure)
  stradl_dat$measure[grep("glob_cc|glob_eff", stradl_dat$dv)] <- "GLOBAL"
  stradl_dat$cohort <- "GS"
  
  comb_dat <- rbind(ukb_dat, stradl_dat)
  comb_dat$abb <- str_wrap(comb_dat$abb, width = 1)
  comb_dat$cohort <- forcats::fct_relevel(comb_dat$cohort, "UKB")
  
  #plot 
  strip_aes <- strip_nested(
    #horizontal strips
    background_x = elem_list_rect(fill = c("midnightblue", "#FFF18D")),
    text_x = elem_list_text(colour = c("white", "black"), face = c("bold", "bold")),
    by_layer_x = T,
    #vertical strips
    background_y = elem_list_rect(fill = "dodgerblue4"),
    text_y = elem_list_text(colour = "white", face = "bold"),
    by_layer_y = T)
  
  f1 <- ggplot(comb_dat, aes(x=abb, y=cohen_d, fill = cohort)) + #either x=abb or x=region
    geom_bar(stat = "identity", width = 1, position = position_dodge(preserve = "total")) +
    geom_errorbar(aes(ymin=cohen_d-sd, ymax=cohen_d+sd),width = .2, position = position_dodge(0.9, preserve = "total")) + #size = 1
    geom_text(aes(label = pval_hfdr_sig, y = 1.1*(cohen_d + sign(cohen_d)*sd)), size = 7, position = position_dodge2(0.9, preserve = "total")) + #hjust=-3,colour="black",,
    labs(y = "Cohen's d",
         title = "Global measures") +
    ylim(-0.3,0.3) +
    scale_fill_manual("Cohort", values = c("#7F867B","#C7C7BB")) + 
    facet_nested(.~measure+abb, scales = "free", strip = strip_aes) + 
    theme(text = element_text(face="bold"),
          axis.text.y = element_text(face="bold", angle = 90),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(), 
          #panel.background = element_rect(size = 1.5), 
          legend.position = "bottom")  

  return(f1)
}

#to plot global alone 
ukb <- ukb_cd[[6]] %>% .[grep("^G", .$measure),] 
stradl <- stradl_cd[[6]] %>% .[.$dv %in% ukb$dv, ]

plt_bestT_cohend <- best_thresh_FUN(ukb, stradl)


#Plot tier cohend ====

#reorganise data just for tiers, for clarity
tier_org_FUN <- function(pval_ls, pval_hfdr_ls){
  dat <- lapply(1:7, function(x){
    tmp <- pval_ls[[x]]
    tmp$thresh <- nm[x]
    tmp$tier <- as.character(tmp$dv) %>% gsub("_local_cc|_local_nef","",.)
    tmp$measure <- as.character(tmp$dv) %>% gsub("t1_|t2_|t3_|t4_","",.)
    tmp$measure <- ifelse(tmp$measure == "local_cc", "LCC", 
                          ifelse(tmp$measure == "local_nef", "LNEF", NA))
    pval_hfdr <- pval_hfdr_ls[[x]]
    tmp1 <- merge(tmp, pval_hfdr[,c("dv","pval_hfdr","pval_hfdr_sig")], by = "dv") 
    
    return(tmp1)})
  names(dat) <- nm
  return(dat)}

ukb_tier <- tier_org_FUN(tier_all_pval_cohend_nodal_cidi_fa, ukb_pval_hfdr) %>% .[[6]] 
stradl_tier <- tier_org_FUN(tier_all_pval_cohend_nodal_cidi_GS_fa, stradl_pval_hfdr) %>% .[[6]] 

plt_tier_FUN <- function(ukb_dat, stradl_dat){
  
  #reorg dat for plotting purposes 
  ukb_dat$cohort <- "UKB"
  stradl_dat$cohort <- "GS"
  ukb_dat$tier <- toupper(ukb_dat$tier)
  stradl_dat$tier <- toupper(stradl_dat$tier)
  
  comb_dat <- rbind(ukb_dat, stradl_dat)
  comb_dat$cohort <- forcats::fct_relevel(comb_dat$cohort, "UKB")
  
  #plot 
  strip_aes <- strip_nested(
    #horizontal strips
    background_x = elem_list_rect(fill = c("midnightblue", "#FFF18D")),
    text_x = elem_list_text(colour = c("white", "black"), face = c("bold", "bold")),
    by_layer_x = T,
    #vertical strips
    background_y = elem_list_rect(fill = "#FFF18D"), #dodgerblue4
    text_y = elem_list_text(colour = "black", face = "bold"), #colour = "white"
    by_layer_y = T)
  
  f1 <- ggplot(comb_dat, aes(x=0, y=cohen_d, fill = cohort)) + #either x=abb or x=region
    geom_bar(stat = "identity", width = 1, position = position_dodge(preserve = "total")) +
    geom_errorbar(aes(ymin=cohen_d-sd, ymax=cohen_d+sd),width = .2, position = position_dodge(0.9, preserve = "total")) + #size = 1
    geom_text(aes(label = pval_hfdr_sig, y = 1.1*(cohen_d + sign(cohen_d)*sd)), size = 7, position = position_dodge2(0.9, preserve = "total")) + #hjust=-3,colour="black",,
    labs(y = "Cohen's d", 
         title = "Tier measures") +
    # ylim(-0.25,0.25) +
    scale_fill_manual("Cohort", values = c("#7F867B","#C7C7BB")) + 
    facet_nested(measure~tier, scales = "free", strip = strip_aes) + 
    theme(text = element_text(face="bold"),
          axis.text.y = element_text(face="bold", angle = 90),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(), 
          #panel.background = element_rect(size = 1.5), 
          legend.position = "bottom")  
  
  return(f1)
}

tier_bestT_cohend <- plt_tier_FUN(ukb_tier, stradl_tier)
tier_bestT_cohend


#Plot nodal tiers in UKB and GS ====

#determine tiers for each subject
#am hardcoding here but can actually convert this section to a function 

#UKB
nodal_dat <- allnodes_cidi_fa[[6]]
deg_dat <- nodal_dat %>% select(ID, contains("deg"))
deg_ls <- split(deg_dat[-1], seq(nrow(deg_dat))) %>% setNames(., deg_dat$ID)

quart_grp <- lapply(deg_ls, function(x) quantile(x)) %>% do.call("rbind", .) %>% rownames_to_column()
colnames(quart_grp) <- c("ID", "t4","t3","t2","t1","max")

deg_quart <- full_join(deg_dat, quart_grp)
quart_dat <- apply(deg_quart[2:86], 2, function(x){ #assign tiers to each node for each subj
  ifelse(x < deg_quart$t3, "t4",
         ifelse(x < deg_quart$t2, "t3",
                ifelse(x < deg_quart$t1, "t2", "t1")))}) %>% 
  as.data.frame() %>% 
  cbind(deg_quart$ID, .)
colnames(quart_dat)[-1] <- gsub("_raw_local_deg","", colnames(quart_dat)[-1])
colnames(quart_dat)[1] <- "ID"

#merge with grouping 
groupings <- readRDS("ukb_mdd_pheno_compiled.rds")
ukb_tier <- merge(groupings[,c(1,6)], quart_dat, , by.x = "connectome_id", by.y = "ID")
colnames(ukb_tier)[1:2] <- c("ID", "group")

#divide into cases and controls and calculate freq
case <- ukb_tier %>% filter(group == 1) %>% select(-group)
case_grp <- apply(case[-1], 2, function(x) return(names(sort(table(x), decreasing = T, na.last = T)[1]))) %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  mutate(., group = 1)
colnames(case_grp) <- c("node","grp", "cidi")

control <-  ukb_tier %>% filter(group == 0) %>% select(-group)
control_grp <- apply(control[-1], 2, function(x) return(names(sort(table(x), decreasing = T, na.last = T)[1]))) %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  mutate(., group = 0)
colnames(control_grp) <- c("node","grp", "cidi")

ukb_grp <- rbind(case_grp, control_grp)
colnames(ukb_grp)[2:3] <- paste0("ukb_", colnames(ukb_grp[2:3]))

#GS 
nodal_dat <- allnodes_cidi_GS_fa[[6]]
deg_dat <- nodal_dat %>% select(ID, contains("deg"))
deg_ls <- split(deg_dat[-1], seq(nrow(deg_dat))) %>% setNames(., deg_dat$ID)

quart_grp <- lapply(deg_ls, function(x) quantile(x)) %>% do.call("rbind", .) %>% rownames_to_column()
colnames(quart_grp) <- c("ID", "t4","t3","t2","t1","max")

deg_quart <- full_join(deg_dat, quart_grp)
quart_dat <- apply(deg_quart[2:86], 2, function(x){ #assign tiers to each node for each subj
  ifelse(x < deg_quart$t3, "t4",
         ifelse(x < deg_quart$t2, "t3",
                ifelse(x < deg_quart$t1, "t2", "t1")))}) %>% 
  as.data.frame() %>% 
  cbind(deg_quart$ID, .)
colnames(quart_dat)[-1] <- gsub("_raw_local_deg","", colnames(quart_dat)[-1])
colnames(quart_dat)[1] <- "ID"

groupings <- readRDS("stradl_mdd_pheno_compiled.rds")
stradl_tier <- merge(groupings[,c(1,3)], quart_dat, by = "ID")
colnames(stradl_tier)[1:2] <- c("ID", "group")

#divide into cases and controls and calculate freq
case <- stradl_tier %>% filter(group == 1) %>% select(-group)
case_grp <- apply(case[-1], 2, function(x) return(names(sort(table(x), decreasing = T, na.last = T)[1]))) %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  mutate(., group = 1)
colnames(case_grp) <- c("node","grp", "cidi")

control <-  stradl_tier %>% filter(group == 0) %>% select(-group)
control_grp <- apply(control[-1], 2, function(x) return(names(sort(table(x), decreasing = T, na.last = T)[1]))) %>% 
  as.data.frame() %>% rownames_to_column() %>% 
  mutate(., group = 0)
colnames(control_grp) <- c("node","grp", "cidi")

stradl_grp <- rbind(case_grp, control_grp)
colnames(stradl_grp)[2:3] <- paste0("stradl_", colnames(stradl_grp[2:3]))

#combine UKB and GS
comb_grp <- cbind(ukb_grp, stradl_grp[-1])

comb_grp$node <- gsub("V","", comb_grp$node)
comb_grp <- comb_grp %>%
  merge(., abb, by.x = "node", by.y = "node_pos", all.x = T) %>% 
  select(-node) %>% 
  select(abb, region, everything())
comb_grp <- comb_grp[order(comb_grp$ukb_grp),]

#melt the cols
tmp <- comb_grp %>% select(-ukb_cidi, -stradl_cidi)
tmp_melt <- pivot_longer(tmp, cols = c("ukb_grp", "stradl_grp"), 
                         names_to = "cohort", names_pattern = "([a-z]{3,6})",
                         values_to = "group")
tmp1 <- comb_grp %>% select(-ukb_grp, -stradl_grp)
tmp1_melt <- pivot_longer(tmp1, cols = c("ukb_cidi", "stradl_cidi"), 
                          names_to = "cohort", names_pattern = "([a-z]{3,6})",
                          values_to = "cidi")
comb_melt <- cbind(tmp_melt, tmp1_melt[,c("cidi")])

comb_melt$group <- toupper(comb_melt$group)
comb_melt$cohort <- as.factor(toupper(comb_melt$cohort))
comb_melt$abb <- as.factor(comb_melt$abb)
comb_melt$cidi <- as.factor(ifelse(comb_melt$cidi == 1, "Cases", "Controls"))

tier_plt <- ggplot(comb_melt, aes(fill = group)) +
  geom_tile(aes(x = fct_inorder(cidi), y = fct_rev(fct_inorder(abb))), color = "white", size = 0.1)+
  labs(title = "Tiers") +
  scale_fill_manual('Tier', values = c("#6D1919","#B30000","#FC8D59","#FDCC8A")) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(.~fct_inorder(cohort)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


#Combine all global and tier plots together ====
#order should be this -> A: plt_bestT_cohend; B: f1; C: tier_bestT_cohend

library(cowplot)
glob_tier <- ggdraw() +
  draw_plot(plt_bestT_cohend, 0, 0.5, 0.5, 0.5, scale = 0.99) + #a - top left
  draw_plot(tier_plt, 0.5, 0, 0.5, 1, scale = 0.99) + #b - long right
  draw_plot(tier_bestT_cohend, 0, 0, 0.5, 0.5, scale = 0.99) + #c - bottom left
  draw_plot_label(c("A", "B", "C"), c(0,0.5, 0), c(1,1, 0.5), size = 20, fontface = "bold")

ggsave(plot = glob_tier, filename = "glob_tier_title.pdf", width = 30, height = 40, units = 'cm', dpi = 300)


