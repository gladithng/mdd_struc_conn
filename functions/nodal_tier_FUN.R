nodal_tier_FUN <- function(nodal_dat, cohort){
  
  deg_dat <- nodal_dat %>% dplyr::select(ID, contains("deg"))
  deg_ls <- split(deg_dat[-1], seq(nrow(deg_dat))) %>% setNames(., deg_dat$ID)
  
  #following Smith et al., each structural connectome was split into four tiers 
  #based on quartiles of the maximum degree. 
  #Tier 1: all nodes with degree greater than 75% of the maximum degree
  #Tier 2: all nodes with degree greater than 50% and up to 75% of the maximum degree
  #Tier 3: all nodes with degree greater than 25% and up to 50% of the maximum degree
  #Tier 4: all nodes with degree greater than 0% and up to 25% of the maximum degree
  
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
  
  #convert df to ls, one ls for each subj
  quart_ls <- split(quart_dat, seq(nrow(quart_dat))) %>% setNames(., quart_dat$ID)
  
  #run function by looping over each subj
  measure <- c("_local_cc", "_local_nef")
  results <- vector("list", 2)
  names(results) <- c("lcc","lnef")
  c <- 1 
  for(i in 1:length(measure)){ #loop over each measure
    subj_id <- nodal_dat$ID
    dat <- nodal_dat[,grep(measure[i], colnames(nodal_dat))]
    dat_ls <- split(dat, seq(nrow(dat))) %>% setNames(., subj_id)
    avg_dat <- data.frame(t1 = NA, t2 = NA, t3 = NA, t4 = NA)
    colnames(avg_dat) <- c(paste0("t1", measure[i]), paste0("t2", measure[i]), 
                           paste0("t3", measure[i]), paste0("t4", measure[i]))
    
    int_FUN <- function(measure_dat, tier_dat){
      for(k in c("t1","t2","t3","t4")){ #loop over each tier
        t_tier <- t(tier_dat[-1]) %>% as.data.frame() %>% rownames_to_column()
        colnames(t_tier) <- c("node_pos","tier")
        grp <- t_tier$node_pos[grep(k, t_tier$tier)] %>% paste0(., collapse = "|") #get tier groupings 
        tmp <- measure_dat[,grep(grp, colnames(measure_dat))] 
        if(length(tmp) == 1){avg_dat[,c] <- tmp} #ie some tiers only got one node
        if(length(tmp) > 1){avg_dat[,c] <- rowMeans(tmp, na.rm = T)}
        c <- c + 1} #average for all nodes in  the tier
      return(avg_dat)}
    
    results[[i]] <- mapply(int_FUN, dat_ls, quart_ls, SIMPLIFY = F, USE.NAMES = T)
  }
  
  lcc <- results[["lcc"]] %>% do.call("rbind", .) %>% as.data.frame() %>% rownames_to_column()
  colnames(lcc)[1] <- "ID"
  
  lnef <- results[["lnef"]] %>% do.call("rbind", .) %>% as.data.frame() %>% rownames_to_column()
  colnames(lnef)[1] <- "ID"
  
  if(cohort == "stradl"){
    comb_results <- merge(lcc, lnef, by = "ID") %>% 
      merge(., groupings[,c(1,3)], by = "ID") %>% 
      merge(., covs, by = "ID")
    colnames(comb_results)[10] <- "group"
    comb_results$group <- ifelse(comb_results$group == 1, "cases", "control")}
  
  if(cohort == "ukb"){
    comb_results <- merge(lcc, lnef, by = "ID") %>% 
      merge(., groupings[,c(1,6)], by.x = "ID", by.y = "connectome_id") %>% 
      merge(., covs, by.x = "ID", by.y = "connectome_id")
    colnames(comb_results)[10] <- "group"
    comb_results$group <- ifelse(comb_results$group == 1, "cases", "control")} #to match cohen d func
  
  return(comb_results)
}
