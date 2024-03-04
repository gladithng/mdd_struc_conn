richclub_w_lists_FUN <- function(mat, reshuffle=c("weights","links"), NR=1000, directed=NULL){ 
  
  #define parameters - just to ensure no confusion in variables
  reshuffle_param <- reshuffle
  NR_param <- NR
  directed_param <- directed
  
  #adapted from tnet package rich club function
  int_func <- function(net_tnet, reshuffle, NR, directed){
    
    #ennsure that the network conform is in the correct format for tnet
    if (is.null(attributes(net_tnet)$tnet)) net <- as.tnet(net_tnet, type = "weighted one-mode tnet") #this will exclude connections that are 0
    if (attributes(net)$tnet != "weighted one-mode tnet") stop("Network not loaded properly")
    
    #to double check if network is directed 
    if(is.null(directed)) {
      tmp <- tnet::symmetrise_w(net, method = "MAX")
      directed <- (nrow(tmp) != nrow(net) | sum(tmp[,"w"]) != sum(net[,"w"]))
      cat(paste("Network is directed: ", directed, sep="")) 
      
      #check if matrix was converted correctly
      tmp <- as.data.frame(as.table(net_tnet)) %>% filter(Freq != 0)
      colnames(tmp) <- c("i","j","w")
      correct <- (nrow(tmp) == nrow(net) | sum(tmp[,"w"]) == sum(net[,"w"]))
      cat(paste("convert correctly: ", correct, sep=""))}
    
    #define prominence and klevels
    prominence <- degree_w(net) #whether you define rich club by degree or strength
    klevels <- max(prominence[,'degree'])
    # klevels <- 50
    
    #internal function to calculate the non-normalised coefficient
    #calculate fraction of weights shared by rich nodes compared w total amount they 
    #could share if they were connected through the strongest links of the network
    `phi` <- function(net){
      
      results <- cbind(k = 1:klevels, rc = NA)
      
      #rank nodes according to weights, with strongest at the top
      wrank <- net[order(net[,"w"],decreasing=T),"w"]
      
      #remove small nodes with degree<k
      for(i in 1:klevels){
        smallnodes <- prominence[prominence[,"degree"]<i,"node"] #return node number
        cutout_net <- net
        cutout_net <- cutout_net[cutout_net[,"i"] %in% smallnodes == F,]
        cutout_net <- cutout_net[cutout_net[,"j"] %in% smallnodes == F,]
        
        #calculate total weight of connections in subset cutout_net
        weights_sum <- sum(cutout_net[,"w"])
        
        #total number of connections in subset cutout_net
        conn_number <- nrow(cutout_net)
        
        #number of connections with max weight in network
        wrank_r <- wrank[1:conn_number]
        
        #weighted rich-club coefficient
        results[i,"rc"] <- weights_sum / sum(wrank_r)
        output <- results[,"rc"]
      }
      return(output)
    }
    
    #calculate the non-normalised coefficient
    ori <- data.frame(k = 1:klevels,rc= phi(net))
    
    #calculate random networks - for each degree k, you have 1000 random networks
    set.seed(1234)
    rand <- matrix(data=0, nrow = nrow(ori), ncol = NR)
    
    for(i in 1:NR) {
      #if(i/10 == round(i/10) )
      #cat(paste("Random network ", i, "/", NR, " @ ", date(), "\n", sep="")) #print this every 10 networks
      rnet <- tnet::rg_reshuffling_w(net, option = reshuffle, directed = F) #this generates a new set of random network each time you run it 
      rand[,i] <- phi(rnet)} #each col represent one network (ncol = # networks)
    
    #define the normalised coefficient
    norm <- data.frame(k = ori[,"k"], rc_norm = 0)
    norm[,"rc_norm"] <- ori[,"rc"]/rowMeans(rand) #take avg across all networks for each klevel
    #norm <- norm[!is.na(norm[,"rc_w_norm"]),]
    
    final <- list(ori,rand,norm)
    names(final) <- c("ori","rand","norm")
    
    return(final)
  }
  
  results <- lapply(mat,int_func, reshuffle=reshuffle_param, NR=NR_param, directed=directed_param)
  
  return(results)
} 
