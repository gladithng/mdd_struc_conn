richclub_permutation_lists_FUN <- function(ls){
  
  #ls: list of list containing both case-control lists w rich club coeff for each threshold
  #m1: compare null distribution and non-normalised observed value
  
  cases_dat <- ls[[1]] 
  control_dat <- ls[[2]]
  N <- dim(cases_dat[["rand"]])[2] #doesn't matter case or control, cuz number of random networks used is the same 
  
  #1. significance of curves for each group - compare to respective null values
  #identify range of k which show rich club organisation (ie when rc_norm is greater than 1 for a continuous range of k)
  cases_k <- as.data.frame(cases_dat[["norm"]]) %>% filter(rc_norm > 1) %>% select(k)
  control_k <- as.data.frame(control_dat[["norm"]]) %>% filter(rc_norm > 1) %>% select(k)
  
  #subset out rc_ori for selected range of k
  cases_ori <- as.data.frame(cases_dat[["ori"]]) %>% .[.$k %in% cases_k$k,] %>% select(rc) 
  control_ori <- as.data.frame(control_dat[["ori"]]) %>% .[.$k %in% cases_k$k,] %>% select(rc) 
  
  #extract the N random networks for selected range of k 
  cases_rand <- as.data.frame(cases_dat[["rand"]]) %>% .[rownames(.) %in% cases_k$k,] %>% 
    split(., seq(nrow(.))) #convert to list of lists, so one list for each k, each k have 1000 elements
  control_rand <- as.data.frame(control_dat[["rand"]]) %>% .[rownames(.) %in% control_k$k,] %>% 
    split(., seq(nrow(.))) #convert to list of lists, so one list for each k 
  
  #% of random null values that exceed observed value for each k
  #observed value should be the non-normalised value
  pval_cases <-  lapply(1:max(seq(cases_rand)), function(x) (sum(cases_rand[[x]] > cases_ori[x,])+1)/(N+1)) #one tailed, not corrected
  eval(parse(text = paste0('pval_cases_m1 <- data.frame(cases_k,pval_m1=unlist(pval_cases))')))
  
  pval_control <-  lapply(1:max(seq(control_rand)), function(x) (sum(control_rand[[x]] > control_ori[x,])+1)/(N+1)) #one tailed, not corrected
  eval(parse(text = paste0('pval_control_m1 <- data.frame(control_k,pval_m1=unlist(pval_control))')))
  
  #2. case-control group difference 
  #generate a null distribution of 10000 random DIFFERENCES for each level of k (ie diff between random null values of group 1 and group 2)
  #pval as the % of null differences that exceeded the observed/original diff for each level of k
  
  #extract common range of k between cases and controls 
  common_k <- 1:nrow(merge(as.data.frame(cases_dat[["norm"]]),as.data.frame(control_dat[["norm"]]),by="k") %>% drop_na()) 
  
  #observed difference in rc_ori between groups for each level of k
  ori_combined <- cbind(as.data.frame(cases_dat[["ori"]])[common_k,c(1:2)], as.data.frame(control_dat[["ori"]])[common_k,c(2)])
  colnames(ori_combined) <- c("k","case_ori","control_ori")
  ori_combined$obsv_diff <- abs(ori_combined$control_ori - ori_combined$case_ori) #calculate observed/original difference
  
  #extract N random null values of group differences for each k (looking at all levels of k here)
  cases_rand <- as.data.frame(cases_dat[["rand"]])[common_k,]  #each row represent k, each col represent one random network 
  control_rand <- as.data.frame(control_dat[["rand"]])[common_k,]
  diff_rand <- abs(control_rand - cases_rand)
  
  #determine significance 
  ls_diff <- diff_rand %>% split(., seq(nrow(.))) #convert to list of lists, so one list for each k
  pval_cases_control <-  lapply(1:max(seq(ls_diff)), function(x) (sum(ls_diff[[x]] > ori_combined$obsv_diff[x])+1)/(N+1)) #one tailed, uncorrected
  eval(parse(text = paste0('pval_cases_control_m1 <- data.frame(k=1:nrow(ori_combined),pval_m1=unlist(pval_cases_control))')))
  
  #extract avg for 1000 random networks - for plotting purposes 
  cases_rand_avg <- as.data.frame(cases_dat[["rand"]]) %>% rowMeans(.) %>% 
    as.data.frame(.) %>% rownames_to_column(.)
  colnames(cases_rand_avg) <- c("k", "cases_rc_rand")
  control_rand_avg <- as.data.frame(control_dat[["rand"]]) %>% rowMeans(.) %>% 
    as.data.frame(.) %>% rownames_to_column(.)
  colnames(control_rand_avg) <- c("k", "control_rc_rand")
  
  #compile results 
  results <-merge(as.data.frame(cases_dat[["ori"]]),as.data.frame(control_dat[["ori"]]),by="k",all=T) %>% 
    merge(.,cases_rand_avg,by="k",all=T) %>% 
    merge(.,control_rand_avg,by="k",all=T) %>% 
    merge(.,as.data.frame(cases_dat[["norm"]]),by="k",all=T) %>% 
    merge(.,as.data.frame(control_dat[["norm"]]),by="k",all=T) %>% 
    merge(.,pval_cases_m1,by="k",all.x=T,all.y=T) %>%
    merge(.,pval_control_m1,by="k",all.x=T,all.y=T) %>% 
    merge(., pval_cases_control_m1,by="k",all.x=T,all.y=T) 
  
  colnames(results) <- c("k","cases_rc_ori","control_rc_ori",
                         "cases_rc_rand","control_rc_rand",
                         "cases_rc_norm","control_rc_norm",
                         "pval_cases_m1","pval_control_m1",
                         "pval_cases_control_m1") 
  
  return(results)
}
