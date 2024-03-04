nodal_measures_indiv_threshprop_FUN <- function(mat, thresh, bin){ 
  
  #"mat" is a list containing 2 lists - cases_mat and control_mat in 3D format (85x85xno_of_subj)
  
  #convert 3d matrix to ls format 
  mat_ls <- lapply(mat, function(x) alply(x,3, .dims = T)) 
  
  #flatten list aka dont put into cases and control 
  mat_ls <- purrr::flatten(mat_ls) #;ts <- head(ts,10)
  
  #prep thresh_FUN
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
  
  #apply mask to each individual network 
  thresh_ls <- lapply(mat_ls, function(x) thresh_FUN(x, thresh))
  
  #change raw mat to igraph object needed to calculate minspantree
  igraph_ls <- lapply(mat_ls, function(x) graph_from_adjacency_matrix(x, weighted = T, mode = "undirected", diag = F))
  
  #calculate minspant tree using raw data and convert back to adjacency
  mst_ls <- lapply(igraph_ls, function(x){
    tmp <- mst(x)
    adj_mst <- as.matrix(as_adjacency_matrix(tmp,attr = "weight"))
    return(adj_mst)})
  
  #add back minspant tree to thresholded matrices and convert to igraph
  add_mst_FUN <- function(mat_dat, mst_dat){
    mat_dat[which(mst_dat != 0)] <- mst_dat[which(mst_dat != 0)]
    return(mat_dat)}
  
  add_mst_ls <- mapply(add_mst_FUN, thresh_ls, mst_ls, SIMPLIFY = F) %>% 
    lapply(., function(x) graph_from_adjacency_matrix(x, weighted = T, mode = "undirected", diag = F))
  
  #check if still disconnected - yes, all returned value of 1
  ts <- lapply(add_mst_ls, function(x) components(x)$no) %>% do.call("rbind", .) %>% as.data.frame() %>% sum(., na.rm = T)
  cat(paste("weighted, component sum: ", ts, sep=""))
  
  #weighted network measures ====
  #adapted from BCT toolbox - for clarity and own understanding of steps
  #measures derived were compared to those derived in matlab to ensure similarity
  
  #to ensure correct matrix used
  if(bin == "no"){ #meaning binary/unweighted
    mat_ls <- add_mst_ls 
    cat("weighted, mat_ls = add_mst_ls") 
    
    #global efficiency 
    glob_eff_FUN <- function(y){
      igraph::E(y)$weight <- igraph::E(y)$weight/max(igraph::E(y)$weight) #need to be scaled between 0 and 1
      igraph::E(y)$weight <- 1/igraph::E(y)$weight #transform weights so that higher weights have shorter distance
      igraph::global_efficiency(y, directed = F)}  #weight automatically used
    glob_eff <- bplapply(mat_ls, glob_eff_FUN)
    df.glob_eff <- glob_eff %>% do.call("rbind",.) %>% as.data.frame(.) %>% rownames_to_column(.) #get in dataframe format
    colnames(df.glob_eff) <- c("ID","glob_eff")
    
    #nodal efficiency
    nodal_eff_FUN <- function(y){
      igraph::E(y)$weight <- igraph::E(y)$weight/max(igraph::E(y)$weight) #need to be scaled between 0 and 1
      igraph::E(y)$weight <- 1/igraph::E(y)$weight #transform weights so that higher weights have shorter distance
      brainGraph::efficiency(y, type = "nodal")}  #weight automatically used
    loc_nef <- bplapply(mat_ls, nodal_eff_FUN)
    df.loc_nef <- loc_nef %>%  do.call("rbind",.) %>% as.data.frame(.) %>% rownames_to_column(.) 
    colnames(df.loc_nef)[1] <- "ID"
    colnames(df.loc_nef)[2:86] <- paste0("V",colnames(df.loc_nef)[2:86])
    colnames(df.loc_nef)[2:86] <- paste0(colnames(df.loc_nef)[2:86],"_raw_local_nef")
    
    #local clustering coefficient
    local_cc_FUN <- function(y){
      igraph::E(y)$weight <- igraph::E(y)$weight/max(igraph::E(y)$weight) #scale mat_ls first
      m <- as.matrix(igraph::as_adjacency_matrix(y, type = "both", attr = "weight"))
      qgraph::clustOnnela(m)[,"clustOnnela"]} 
    loc_cc <- bplapply(mat_ls, local_cc_FUN)
    df.loc_cc <- loc_cc %>%  do.call("rbind",.) %>% as.data.frame(.) %>% rownames_to_column(.) 
    colnames(df.loc_cc)[1] <- "ID"
    colnames(df.loc_cc)[2:86] <- paste0(colnames(df.loc_cc)[2:86],"_raw_local_cc")
    
    #global clustering coefficient
    df.glob_cc <- data.frame(ID = df.loc_cc[1], glob_cc = rowMeans(df.loc_cc[-1], na.rm = T))
    glob_cc <- setNames(split(df.glob_cc, seq(nrow(df.glob_cc))), df.glob_cc$ID)
    
    #local degree - needed for tier-based calculations
    local_deg_FUN <- function(y){igraph::degree(y)}
    loc_deg <- bplapply(mat_ls, local_deg_FUN)
    df.loc_deg <- loc_deg %>%  do.call("rbind",.) %>% as.data.frame(.) %>% rownames_to_column(.) 
    colnames(df.loc_deg)[1] <- "ID"
    colnames(df.loc_deg)[2:86] <- paste0("V",colnames(df.loc_deg)[2:86])
    colnames(df.loc_deg)[2:86] <- paste0(colnames(df.loc_deg)[2:86],"_raw_local_deg")
    
    #compile nodal measures - both list and df format
    obs_results <- list(setNames(list(glob_cc, glob_eff, loc_cc, loc_nef, loc_deg), c("glob_cc","glob_eff", "loc_cc", "loc_nef", "loc_deg")),
                        setNames(list(df.glob_cc, df.glob_eff, df.loc_cc, df.loc_nef, df.loc_deg), c("df.glob_cc", "df.glob_eff", "df.loc_cc", "df.loc_nef", "df.loc_deg")))
    names(obs_results) <- c("ls","df")
  }
  
  return(obs_results)
}
