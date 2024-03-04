nodal_cohend_wperm_FUN <- function(dat, thresh, cohort, ls_covs, perm_N){
  
  dat$group <- ifelse(dat$group == "cases", 1, 0)
  dat$group <- as.factor(dat$group)
  dat$group <- relevel(dat$group, ref = "0") #use control as reference level
  
  #1 - run regression model
  
  #define regression function to calculate std_beta, cohens_d, tval, empirical p
  reg_FUN <- function(reg_dat){
    
    ls.dep <- colnames(reg_dat)[grep("^glob_eff$|^glob_cc$|_local_cc$|_local_nef$",colnames(reg_dat))]
    
    ls.factor <- colnames(reg_dat)[grep("group", colnames(reg_dat))]
    ls.dep.factor <- expand.grid(ls.dep, ls.factor, stringsAsFactors = F)
    ls.models <- data.frame(dependent = ls.dep.factor$Var1, 
                            factor = ls.dep.factor$Var2, 
                            covs = "", 
                            stringsAsFactors = F)
    ls.models$covs <- ls_covs
    
    int_FUN <- function(ls, int_dat){
      
      depvar <- as.character(ls[1])
      indpvar <- as.character(ls[2])
      cov <- as.character(ls[3])
      
      #continuous IV, categorical DV, continuous and categorical covariates 
      eq <- paste0('scale(',depvar,')~', cov,'+', indpvar) 
      fit <- lm(as.formula(eq),dat = int_dat, na.action = na.exclude)
      tb <- tidy(fit) 
      
      #calculate cohen's d - Nakagawa & Cuthill, 2007
      tstat <- tb[[nrow(tb),4]] #tstat for group var - controls is used as the reference group 
      n1 <- table(int_dat$group)[[1]] #sample size for group 1 (control)
      n2 <- table(int_dat$group)[[2]] #sample size for group 2 (cases)
      df <- summary(fit)$fstatistic[[3]] #df for whole model
      cohens_d <- (tstat*(n1+n2))/(sqrt(n1*n2)*sqrt(df))
      
      #save model output
      stats <- tail(tb,1)
      mod_result <-  data.frame(dv = depvar, iv = indpvar, stats[,c(2:5)], cohen_d = cohens_d)
      colnames(mod_result) <- c("dv", "iv", "beta", "sd", "tval", "pval_mod", "cohen_d")
      
      return(mod_result) 
    } #end of int_FUN
    
    reg_result <- pbapply(ls.models, 1, int_FUN, reg_dat) %>% do.call("rbind", .)
    return(reg_result)
  } #end of reg_FUN
  
  #run reg_FUN for original data
  ori_reg <- reg_FUN(dat)
  ori_tval <- abs(ori_reg$tval) 
  
  #2 - do permutation to determine final p-value "pval_perm" 
  
  #define permutation function
  #to generate tvals aka run reg analysis for 1000 permutations of original data
  perm_FUN <- function(n, ori_dat){
    
    perm_group <- sample(ori_dat$group)
    perm_dat <- cbind(ori_dat, perm_group) %>% select(-group) #so you're only left with perm_group #select(ID, group, perm_group)
    colnames(perm_dat)[grep("perm_group", colnames(perm_dat))] <- "group" #rename to make things easier
    
    #run reg_FUN for perm_dat 
    perm_reg <- reg_FUN(perm_dat)
    perm_tval <- abs(perm_reg$tval)
    
    return(perm_tval)
  }
  
  #run perm_FUN perm_N times 
  perm_vals <- foreach(i = 1:perm_N) %dorng% {perm_FUN(i, dat)} #end up with a list
  perm_vals <- perm_vals %>% do.call("rbind", .) %>% t(.) #each row represnt a measure, each col represent a permutation
  
  #determine pval by counting number of tstat that exceed original tstat for each measure 
  pdat <- data.frame(type = colnames(dat)[grep("^glob_eff$|^glob_cc$|_local_cc$|_local_nef$",colnames(dat))], pval_perm = NA)
  j <- 1
  for(i in 1:nrow(perm_vals)){
    pdat[j,2] <- (sum(perm_vals[i,] > ori_tval[i])+1)/(perm_N+1)
    j <- j + 1}
  
  #combine pval_perm with ori_reg 
  comb_dat <- merge(ori_reg, pdat, by.x = "dv", by.y = "type")
  
  #3 - rename nodes for ease of interpretation
  #"nodes_85" and "abb" variables must be in the global env
  nodes <- gsub("(V)([0-9]{1,2})(\\_.+)", "\\2", comb_dat$dv[grep("local", comb_dat$dv)]) %>% 
    unique() %>% as.data.frame()
  region <- merge(nodes_85, nodes, by.x = "node_pos", by.y=".") 
  region$node_pos <- as.numeric(region$node_pos)
  region <- region[order(region$node_pos),]
  
  comb_dat$region <- comb_dat$dv
  comb_dat <- comb_dat %>% 
    tidyr::extract(., col = "region", into = c(NA,"node_pos",NA,NA), "(V)([0-9]{1,2})(\\_)([a-z]{3,4}\\_[a-z]{5}\\_[a-z]{2,3})") %>% 
    merge(., region, by = "node_pos", all.x = T) %>% 
    merge(., abb, by = c("node_pos", "region"), all.x = T) %>% 
    select(-node_pos) %>% 
    select(abb, region, everything())
  
  return(comb_dat) 
}
