## PLOT GGSEG AND P-PLOT FOR NODAL MEASURES ##
#for best performing threshold only aka 35% 

rm(list=ls())
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr);library(plyr)
library(dplyr);library(ggplot2);library(magrittr);library(ggpubr);library(cowplot)
library(ggrepel);library(tidytext);library(ggh4x);library(ggseg);
library(viridis);library(patchwork);library(gridExtra)

#Load data ====

nodes_85 <- read.csv("fslabels_85node.txt", header=F) %>% rownames_to_column(.)
colnames(nodes_85) <- c("node_pos","region")

atlas_abb <- readxl::read_excel("nodes_85_numbered.xlsx", sheet = "atlas")
abb <- readxl::read_excel("C:/Users/gladi/OneDrive/Desktop/connectomes/analysis/stats/nodes_85_numbered.xlsx")

#empirical pvals for global, tier and nodal 
load("ukb_faonly_cidi_PT_nodal_pval.RData")
load("stradl_faonly_cidi_PT_nodal_pval.RData")

#hfdr-corrected pvals for global, tier and nodal
ukb_pval_hfdr <- readRDS("ukb_nodal_tier_hfdr_pval.rds")
stradl_pval_hfdr <- readRDS("stradl_nodal_tier_hfdr_pval.rds")

#Merge data with pval_hfdr ====
nm <- c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40")
merge_FUN <- function(pval_ls, pval_hfdr_ls){
  dat <- lapply(1:7, function(x){
    tmp <- pval_ls[[x]]
    pval_hfdr <- pval_hfdr_ls[[x]]
    tmp1 <- merge(tmp, pval_hfdr[,c("dv","pval_hfdr","pval_hfdr_sig")], by = "dv") #merge with pval_hfdr 
    return(tmp1)})
  names(dat) <- nm
  return(dat)}

ukb_df <- merge_FUN(all_pval_cohend_nodal_cidi_fa, ukb_pval_hfdr) %>% .[[6]]
stradl_df <- merge_FUN(all_pval_cohend_nodal_cidi_GS_fa, stradl_pval_hfdr) %>% .[[6]]


#Plot cohen's d using ggseg ====

org_FUN <- function(nodal_df, sig = T, dk = T, aseg = F){
  nodal_df$hemi <- NA
  nodal_df$hemi[grep("^L", nodal_df$abb)] <- "left"
  nodal_df$hemi[grep("^R", nodal_df$abb)] <- "right"
  nodal_df$hemi[grep("bSTEM", nodal_df$abb)] <- "midline" #following ggseg aseg atlas: aseg$data
  
  nodal_df$atlas <- "aseg"
  nodal_df$atlas[grep("^ctx", nodal_df$region)] <- "dk" #desikian
  
  #format names to suit this 
  #https://github.com/ggseg/ggseg
  nodal_df <- merge(atlas_abb[c(2,4)], nodal_df, by = "region")
  nodal_df$atlas_region <- gsub("Left-","", nodal_df$atlas_region)
  nodal_df$atlas_region <- gsub("Right-","", nodal_df$atlas_region)
  nodal_df$atlas_region <- gsub("ctx-lh-","", nodal_df$atlas_region)
  nodal_df$atlas_region <- gsub("ctx-rh-","", nodal_df$atlas_region)
  
  #if plotting sig regions only based on pval_perm
  if(sig == T){nodal_df <- nodal_df %>% filter(pval_perm < 0.05)}
  
  #to separate into dk and aseg
  if(dk == T){
    dat <- nodal_df %>% 
      filter(atlas == "dk") %>% 
      select(atlas_region, dv, hemi, cohen_d) 
    colnames(dat)[1] <- "region"}
  
  if(aseg == T){
    dat <- nodal_df %>% 
      filter(atlas == "aseg") %>% 
      select(atlas_region, dv, hemi, cohen_d)
    colnames(dat)[1] <- "region"}
  
  return(dat)
}

#Prep dat for cortical regions (dk)
ukb_dk_dat <- org_FUN(ukb_df, sig = F, dk = T, aseg = F)
stradl_dk_dat <- org_FUN(stradl_df, sig = F, dk = T, aseg = F)

#Prep dat for subcortical regions (aseg)
#aseg atlas have no accumbens area ...
ukb_aseg_dat <- org_FUN(ukb_df, sig = F, dk = F, aseg = T)
stradl_aseg_dat <- org_FUN(stradl_df, sig = F, dk = F, aseg = T)

#Plot for each type of nodal measure (LCC and LNEF)
plot_ggseg_FUN <- function(dk_dat, aseg_dat, cohort, legend){
  
  #cortical, LCC
  f1_dat <- dk_dat[grep("_local_cc$", dk_dat$dv), -grep("dv", colnames(dk_dat))]
  f1 <- f1_dat %>% ggplot() +
    geom_brain(mapping = aes(fill = as.numeric(cohen_d)), atlas = dk, 
               position = position_brain(hemi~side)) + #position = position_brain(hemi~side); .~side+hemi
    scale_fill_gradient2(low = "#AE1D24", mid = "white", high ="#0571B1", midpoint = 0,
                         limits = c(-0.2, 0.2), breaks = c(-0.2, 0, 0.2)) +
  labs(title = "CC",
       fill = "Cohen's d") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          legend.position = "none") 
  
  #cortical, LNEF
  f2_dat <- dk_dat[grep("_local_nef$", dk_dat$dv), -grep("dv", colnames(dk_dat))]
  f2 <- f2_dat %>% ggplot() +
    geom_brain(mapping = aes(fill = as.numeric(cohen_d)), atlas = dk, 
               position = position_brain(hemi~side)) + #position = position_brain(hemi~side); .~side+hemi
    scale_fill_gradient2(low = "#AE1D24", mid = "white", high ="#0571B1", midpoint = 0, 
                         limits = c(-0.2, 0.2), breaks = c(-0.2, 0, 0.2)) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) + #https://ggplot2.tidyverse.org/reference/guide_colourbar.html
    labs(title = "NEFF", fill = "Cohen's d") +
    theme_minimal() +
    theme(axis.text = element_blank(), 
          legend.title = element_text(face = "bold"),
          legend.position = "bottom") 
  
  #subcortical, LCC
  f3_dat <- aseg_dat[grep("_local_cc$", aseg_dat$dv), -grep("dv", colnames(aseg_dat))]
  f3 <- f3_dat %>% ggplot() +
    geom_brain(mapping = aes(fill = as.numeric(cohen_d)), atlas = aseg, side = "coronal") +
    scale_fill_gradient2(low = "#AE1D24", mid = "white", high ="#0571B1", midpoint = 0,
                         limits = c(-0.2, 0.2), breaks = c(-0.2, 0, 0.2)) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          legend.position = "none",
          legend.title = element_text(face = "bold"))
  
  #subcortical, LNEF
  f4_dat <- aseg_dat[grep("_local_nef$", aseg_dat$dv), -grep("dv", colnames(aseg_dat))]
  f4 <- f4_dat %>% ggplot() +
    geom_brain(mapping = aes(fill = as.numeric(cohen_d)), atlas = aseg, side = "coronal") +
    scale_fill_gradient2(low = "#AE1D24", mid = "white", high ="#0571B1", midpoint = 0, 
                         limits = c(-0.2, 0.2), breaks = c(-0.2, 0, 0.2)) +
    theme_minimal() +
    theme(axis.text = element_blank(),
          # text = element_text(face = "bold"), 
          legend.position = "none", #bottom
          legend.title = element_text(face = "bold"))
  
  #Combine all plots
  comb <- (f1|f3)/(f2|f4) + 
    patchwork::plot_annotation(title = cohort, 
                               theme = theme(plot.title = element_text(face = "bold", hjust = 0.5),
                                             legend.position =  legend)) +
    patchwork::plot_layout(guides = "collect")
  return(comb)
}

ukb_ggseg <- plot_ggseg_FUN(ukb_dk_dat, ukb_aseg_dat, cohort = "UKB", legend = "none")
stradl_ggseg <- plot_ggseg_FUN(stradl_dk_dat, stradl_aseg_dat, cohort = "GS", legend = "none")

#combine ggseg
# comb_ggseg <- ggpubr::ggarrange(ukb_ggseg, stradl_ggseg, ncol=2, legend = "none") %>% 
#   gridExtra::grid.arrange(get_legend(ukb_ggseg), heights = unit(c(70, 10), "mm"),
#                           top = grid::textGrob(expression(bold("Nodal measures")), x = 0.05, y = 6, hjust = 0)) 

#Plot correlation between UKB and GS for all the nodal measures ====

ukb_loc <- ukb_df[grep("_local_cc$|_nef$", ukb_df$dv), c("dv", "cohen_d")]
colnames(ukb_loc)[2] <- "ukb_cohen"
stradl_loc <- stradl_df[grep("_local_cc$|_nef$", stradl_df$dv), c("dv", "cohen_d")]
colnames(stradl_loc)[2] <- "stradl_cohen"

comb_loc <- merge(ukb_loc, stradl_loc, by = "dv")
comb_loc$measure[grepl("nef", comb_loc$dv)] <- "NEFF"
comb_loc$measure[grepl("cc", comb_loc$dv)] <- "CC"

loc_corr <- ggplot(comb_loc, aes(x = ukb_cohen, y = stradl_cohen)) + 
  geom_point(aes(color = measure)) + 
  geom_smooth(method = lm, se = F, color = "black") +
  stat_cor(method = "pearson", label.x.npc = 0, label.y.npc = "top") +
  scale_color_brewer("Measure", palette = "Dark2") +
  labs(title = "Correlation between UKB and GS",
       x = "UKB Cohen's d", 
       y = "GS Cohen's d") + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), 
        legend.position = "bottom") #r=0.3 and sig

#Plot pval thresholds for nodal measures ====
#adapted from Shen's code 

p_plot<-function(TargetResult, color_theme = 'custom', shape_sig = T,
                 outputpath = NA, category.input, labels_annot, plot_title = NA,
                 add_category_name = F, fig_size = c(24,20), y_lim = 2.5){
  
  #change to correct class 
  TargetResult <- TargetResult %>% 
    select(-grp) %>% 
    mutate_at(vars(beta:pval_perm, pval_hfdr), as.numeric)
  
  #add label for significant results to have a different shape
  TargetResult$shape <- 1 #empty circle
  TargetResult$sig <- ""
  if (shape_sig == T){
    cols <- grepl("local",TargetResult$dv) & TargetResult$pval_hfdr < 0.05
    TargetResult$shape[cols] <- 19 #solid circle
    TargetResult$sig[cols] <- '*'
  }else{}
  
  #set color theme
  if(length(color_theme) > 1){
    cl.theme <- color_theme
  }else if (color_theme == 'custom'){
    cl.theme <- c("#6D1919","#B30000","#FC8D59","#FDCC8A")
  }else{}
  
  #identify significance line
  hfdr_pval <- -log10(0.05)
  hfdr_col <- "#0000CC"
    
  #reorder result table by category and add another column of category
  tmp.result.table <- TargetResult
  tmp.result.table$category =''
  add_category <- function(ls.capture,category.name,targetMat,col.tocap){
    loc.tocap <- grep(ls.capture,targetMat[,col.tocap])
    targetMat[loc.tocap,'category'] <- category.name
    return(targetMat)}
  for(i in 1:nrow(category.input)){
    if(i==1){
      tmp.result.table <- add_category(category.input[i,1], category.input[i,2],
                                       targetMat = TargetResult, col.tocap = 'dv')
    }else{
      tmp.result.table <- add_category(category.input[i,1], category.input[i,2],
                                       targetMat = tmp.result.table, col.tocap = 'dv')
    }
  }
  TargetResult <- tmp.result.table
  
  #add labels
  if(sum(!is.na(labels_annot)) > 0){
    TargetResult$labels <- rep(99999,nrow(TargetResult))
    for(i in 1:nrow(labels_annot)){
      TargetResult$labels[grep(labels_annot[i,1],TargetResult$dv)] <- labels_annot[i,2]}
    TargetResult$labels[TargetResult$labels==99999] <- ''
  }else{
    TargetResult$labels <- ''
  }
  
  #add a column to fix the order of phenotypes
  TargetResult$ord <- 1:nrow(TargetResult)
  
  #set x axis to print categories
  axis.set <- TargetResult %>% 
    dplyr::group_by(category) %>% 
    dplyr::summarize(center = (max(ord) + min(ord)) / 2)
  
  #make the plot
  fig <- ggplot(TargetResult, aes(x = ord, y = -log10(pval_hfdr), label = labels)) +
    geom_point(size = 2, shape = TargetResult$shape, stroke = 2,
               aes(colour = category)) +
    scale_colour_manual(name = "Tier", values = cl.theme)+
    geom_text_repel(box.padding = unit(1,'lines'), segment.size = 0.2,
                    max.iter = 1000, ylim = c(0, Inf)) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face='bold'),
          axis.text.y = element_text(face = "bold"),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "grey"),
          axis.line.x = element_blank(),
          plot.title = element_text(lineheight = 1,face='bold'), 
          legend.position = "top",
          legend.title = element_text(face = "bold")) +
    geom_hline(yintercept = 0 , color = "grey", size=0.5)+
    geom_hline(yintercept = hfdr_pval, linetype="dashed", color = hfdr_col, size = 1) + #local fdr pval
    annotate("text", x = -Inf, y = hfdr_pval+0.05, label = expression(bolditalic('pFDR')), color = hfdr_col, hjust = 0, size = 3) +
    ylab('-log10(p)')+
    ylim(c(0, y_lim))
  
  if(!is.na(plot_title)){
    fig <- fig + ggtitle(plot_title)}
  if(add_category_name==T){
    fig <- fig + scale_x_continuous(label = axis.set$category, 
                                    breaks = axis.set$center)
  }else{
    fig=fig+theme(axis.text.x=element_blank())
  }
  
  #print figure to a file
  if(!is.na(outputpath)){
    ggsave(plot = fig,filename = outputpath,
           width = fig_size[1], height = fig_size[2],
           units = 'cm', dpi = 300)
  }else{}
  
  return(fig)
}

#Based on UKB (since no sig results for GS)

#category
ls.cat <- ukb_pval_hfdr %>% .[[6]] %>% 
  .[grep("_raw_local_nef$", .$dv),] %>% 
  dplyr::select(dv, grp)
ls.cat$grp <- toupper(ls.cat$grp)
colnames(ls.cat) <- c("dependent", "cat")

#labels
ls.lab <- ukb_df$dv[grep("_local_nef$", ukb_df$dv)] %>% as.data.frame() #_local_cc$
colnames(ls.lab) <- "dependent" 

#rename nodes 
nodes <- gsub("(V)([0-9]{1,2})(\\_.+)", "\\2", ls.lab$dependent[grep("local",ls.lab$dependent)]) %>% 
  unique() %>% as.data.frame()
region <- merge(nodes_85, nodes, by.x = "node_pos", by.y=".") 
region$node_pos <- as.numeric(region$node_pos)
region <- region[order(region$node_pos),]

ls.lab$region <- ls.lab$dependent
ls.lab <- ls.lab%>% 
  tidyr::extract(., col = "region", into = c(NA,"node_pos",NA,NA), "(V)([0-9]{1,2})(\\_)([a-z]{3,4}\\_[a-z]{5}\\_[a-z]{2,3})") %>% 
  merge(., region, by = "node_pos", all.x = T) %>% 
  merge(., abb, by = c("node_pos", "region"), all.x = T) %>% 
  select(-node_pos, -region, dependent, abb) 

#label only sig regions
tmp <- ukb_df$dv[-grep("\\*", as.character(ukb_df$pval_hfdr_sig))]
ls.lab$abb[ls.lab$dependent %in% tmp] <- ""

#Run function
ukb_plt_dat <- ukb_df[grep("_local_nef$", ukb_df$dv),]
ukb_pval <- p_plot(ukb_plt_dat, color_theme = "custom", shape_sig = T, 
                   category.input = ls.cat, labels_annot = ls.lab, outputpath = NA, 
                   plot_title = "UKB pFDR for NEFF", add_category_name = T)


#Combine ggseg plots and corrplot ====
#https://wilkelab.org/cowplot/articles/shared_legends.html !!!!!
legend_dat <- ukb_dk_dat[grep("_local_cc$", ukb_dk_dat$dv), -grep("dv", colnames(ukb_dk_dat))]
legend_plt <- legend_dat %>% ggplot() +
  geom_brain(mapping = aes(fill = as.numeric(cohen_d)), atlas = dk, 
             position = position_brain(hemi~side)) + #position = position_brain(hemi~side); .~side+hemi
  scale_fill_gradient2(low = "#AE1D24", mid = "white", high ="#0571B1", midpoint = 0,
                       limits = c(-0.2, 0.2), breaks = c(-0.2, 0, 0.2)) +
  labs(fill = "Cohen's d") +
  theme(legend.position = "bottom", 
        legend.title = element_text(face = "bold")) 

legend_b <- get_legend(legend_plt + 
                         guides(color = guide_legend(nrow = 1)) +
                         theme(legend.position = "bottom"))

p1 <- plot_grid(ukb_ggseg, stradl_ggseg, nrow = 1)
p2 <- plot_grid(p1, legend_b, ncol = 1, rel_heights = c(1, .2))
p3 <- plot_grid(loc_corr, ukb_pval, ncol = 2, labels = c("B", "C"))
p4 <- plot_grid(p2, p3, ncol = 1, labels = c("A", ""))
