hfdr_FUN <- function(nodal_df, nodal_pval, tier_pval){ 
  #nodal_df is allnodes_cidi_fa; tier_pval is tier_all_pval_cohend_nodal_cidi_fa; nodal_pval is all_pval_cohend_nodal_cidi_fa
  
  #to determine which nodes fall into what tiers
  nodal_dat <- nodal_df
  deg_dat <- nodal_dat %>% dplyr::select(ID, contains("deg"))
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
  colnames(quart_dat)[-1] <- gsub("raw_local_deg","", colnames(quart_dat)[-1])
  colnames(quart_dat)[1] <- "ID"
  
  ts <- quart_dat
  
  #to assign nodes into their tiers
  # names(which.max(table(ts$V1_)))
  grp <- apply(ts[-1], 2, function(x) return(names(sort(table(x), decreasing = T, na.last = T)[1]))) %>% 
    as.data.frame() %>% rownames_to_column()
  colnames(grp) <- c("node","grp")
  freq <- table(grp$grp) %>% as.data.frame()
  
  #to prep unadj pvals
  nodal <- nodal_pval
  tier <- tier_pval
  pval <- rbind(nodal, tier) %>% dplyr::select(dv, pval_perm)
  pval$dv <- as.character(pval$dv)
  
  pval$grp <- ""
  pval$grp[grep(paste0(grp$node[grep("t1", grp$grp)], collapse = "|"), pval$dv)] <- "t1"
  pval$grp[grep(paste0(grp$node[grep("t2", grp$grp)], collapse = "|"), pval$dv)] <- "t2"
  pval$grp[grep(paste0(grp$node[grep("t3", grp$grp)], collapse = "|"), pval$dv)] <- "t3"
  pval$grp[grep(paste0(grp$node[grep("t4", grp$grp)], collapse = "|"), pval$dv)] <- "t4"
  
  pval_tb <- pval
  pval_tb$pval_hfdr <- as.numeric(99999)
  
  #correct glob level
  pval_tb$pval_hfdr[grep("^glob_cc$|^glob_eff$", pval_tb$dv)] <- p.adjust(pval_tb$pval_perm[grep("^glob_cc$|^glob_eff$", pval_tb$dv)], method = "fdr")
  str <- paste0(pval_tb$dv[pval_tb$pval_hfdr < 0.05], collapse = "|") #select cols that <0.05 for glob correction 
  suffix_keep <- stringr::str_split_i(str, "_", i = 2) %>% unique()
  
  #correct tier level
  for(i in suffix_keep){
    if(i == "eff"){kp <- "^t\\d+_local_nef$"}
    else{kp <- "^t\\d+_local_cc$"}
    pval_tb$pval_hfdr[grep(kp, pval_tb$dv)] <- p.adjust(pval_tb$pval_perm[grep(kp, pval_tb$dv)], method = "fdr")
  }
  str <- pval_tb$dv[pval_tb$pval_hfdr < 0.05] %>% .[grep("^t", .)] #select cols that <0.05 for tier correction 
  suffix_keep <- stringr::str_split_i(str, "_", i = 3) %>% unique()
  
  #correct nodal level - group acc to nodes in each tier 
  for(i in suffix_keep){
    if(i == "nef"){kp <- "_raw_local_nef$"}
    else{kp <- "_raw_local_cc$"}
    
    for(i in c("t1", "t2", "t3", "t4")){
      tmp <- grepl(kp, pval_tb$dv) & grepl(i, pval_tb$grp)
      pval_tb$pval_hfdr[tmp] <- p.adjust(pval_tb$pval_perm[tmp], method = "fdr")
    }
  }
  
  #insert notation to indicate significance
  pval_tb$pval_hfdr_sig <- ""
  pval_tb$pval_hfdr_sig <- ifelse(pval_tb$pval_hfdr < 0.05, "*", pval_tb$pval_hfdr_sig)
  
  return(pval_tb)
}
