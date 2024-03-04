## PLOT CURVES FOR RICH CLUB ##

#used plotly to allow for the plotting of two y-axes

rm(list=ls())
options(scipen=999)
library(tidyverse); library(broom);library(dplyr);library(stringr)
library(ggplot2);library(reshape);library(reshape2)
library(ggpubr);library(ggh4x);library(cowplot);library(plotly)

#load data 
load("ukb_faonly_cidi_richcoefs_pval.RData")

#organise data suitable for plotly
org_dat_FUN <- function(grp_dat){ 
  
  nm <- sapply(grp_dat,`[[`,1) %>% lapply(., nrow)
  
  dat <- do.call("rbind",sapply(grp_dat,`[[`,1)) %>%
    select(k, contains(c("_norm","pval"))) %>%
    mutate(., thresh = rep(c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"), nm)) 
  # filter(k %in% 1:50)
  colnames(dat)[6] <- "pval_casescontrol_m1" #to facilitate melt later
  
  dat_melt <- dat %>%
    pivot_longer(., cols = 2:3, names_to = c("group","remove"),
                 names_pattern = "([a-z]{5,7})(\\_rc\\_.+)",
                 values_to = "rc_norm") %>%
    mutate(., pval = ifelse(.$group == "cases", .$pval_cases_m1, .$pval_control_m1)) %>%
    select(-remove)
  
  # dat_melt <- cbind(tmp_ori, tmp_rand, tmp_norm)
  # dat_melt <- dat_melt[, !duplicated(colnames(dat_melt))]
  
  # dat_melt <- dat %>%
  # reshape(., idvar = c("k","thresh"), varying = 2:3, v.names = "rc_norm",
  #         direction = "long", timevar = "group", times = c("cases","control")) %>%
  # mutate(., pval = c(.[c(1:150), "pval_cases_m1"], .[c(1:150), "pval_control_m1"]))
  
  #to determine coordinates for shaded region ===
  dat_melt$xmin <- ""
  dat_melt$xmax <- ""
  dat_melt$yrect <- ""
  dat_melt$ylimit <- ""
  result <- dat_melt[0,]
  
  for(i in c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40")){
    tmp <- dat_melt[grepl(i, dat_melt$thresh),] #!is.na(dat_melt$pval) & 
    cases_idx <- grepl("cases", tmp$group)
    control_idx <- grepl("control", tmp$group)
    
    #rich club is defined as rc_norm > 1 (where does the significance come in then..?)
    df <- tmp[cases_idx,] %>% filter(.$rc_norm > 1 & .$pval_cases_m1 < 0.05)
    if(nrow(df) == 0){cases_min <- 0; cases_max <- 0; cases_yrect <- 0}
    if(nrow(df) == 1){cases_min <- head(df$k, 1); cases_max <- tail(df$k, 1); cases_yrect <- 1.05*max(df$rc_norm)}
    if(nrow(df) > 1){
      cases_yrect <- 1.05*max(df$rc_norm) #this only takes into account the max rc_norm that >1 and <0.05
      #test for continuous range
      cont_cases <- split(df$k, cumsum(c(1, diff(df$k) != 1)))
      lg_cases <- lapply(cont_cases, length) >= 3 #arbitrary number - continuous for at least 3 degrees
      xvals_cases <- cont_cases[lg_cases] %>% unlist() %>% unname() #select logical T and unlist the list
      if(is.null(xvals_cases)){cases_min <- 0; cases_max <- 0}
      if(!is.null(xvals_cases)){cases_min <- xvals_cases[1];cases_max <- xvals_cases[length(xvals_cases)]}}
    
    df <- tmp[control_idx,] %>% filter(.$rc_norm > 1 & .$pval_control_m1 < 0.05)
    if(nrow(df) == 0){control_min <- 0; control_max <- 0; control_yrect <- 0}
    if(nrow(df) == 1){control_min <- head(df$k, 1); control_max <- tail(df$k, 1); control_yrect <- 1.05*max(df$rc_norm)}
    if(nrow(df) > 1){
      control_yrect <- 1.05*max(df$rc_norm) #this only takes into account the max rc_norm that >1 and <0.05
      #test for continuous range
      cont_control <- split(df$k, cumsum(c(1, diff(df$k) != 1)))
      lg_control <- lapply(cont_control, length) >= 3 #arbitrary number - continuous for at least 3 degrees
      xvals_control <- cont_control[lg_control] %>% unlist() %>% unname() #select logical T and unlist the list
      if(is.null(xvals_control)){control_min <- 0; control_max <- 0}
      if(!is.null(xvals_control)){control_min <- xvals_control[1]; control_max <- xvals_control[length(xvals_control)]}}
    
    #select higher of the 2 values as ylim to ensure standardisation
    df <- tmp[cases_idx,]; cases_ylimit <- max(df$rc_norm, na.rm = T)
    df <- tmp[control_idx,]; control_ylimit <- max(df$rc_norm, na.rm = T)
    lim <- max(cases_ylimit, control_ylimit, na.rm = T)
    
    lim_rect <- min(cases_yrect, control_yrect, na.rm = T)
    
    tmp$xmin <- ifelse(tmp$group == "cases", cases_min, control_min)
    tmp$xmax <- ifelse(tmp$group == "cases", cases_max, control_max)
    # tmp$yrect <- ifelse(tmp$group == "cases", cases_yrect, control_yrect)
    # tmp$ylimit <- ifelse(tmp$group == "cases", cases_ylimit, control_ylimit)
    tmp$yrect <- ifelse(tmp$group == "cases", lim_rect, lim_rect)
    tmp$ylimit <- ifelse(tmp$group == "cases", lim_rect, lim_rect)
    
    result <- rbind(result, tmp)}
  
  dat_melt <- result
  
  #to determine significance for between group differences  ===
  # dat_melt$mdd <- mdd
  dat_melt$both_sig <- ifelse(dat_melt$pval_cases_m1 < 0.05 & dat_melt$pval_control_m1 < 0.05, 1,NA)
  result1 <- dat_melt[0,]
  
  for(i in c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40")){
    tmp <- dat_melt[grepl(i, dat_melt$thresh),]
    idx_k <- tmp$k[!is.na(tmp$pval_cases_m1) & !is.na(tmp$pval_control_m1)]
    tmp$grp_sig[tmp$k %in% idx_k] <- ifelse(tmp$pval_casescontrol_m1[tmp$k %in% idx_k] < 0.05, "*", "")
    result1 <- rbind(result1, tmp)}
  
  dat_melt <- result1
  
  #merge with ori and rand data - for plotting purposes if using sec.axis ===
  for(i in c("ori", "rand")){
    tmp <- do.call("rbind",sapply(grp_dat,`[[`,1)) %>%
      select(k, contains(c(i,"pval"))) %>%
      mutate(., thresh = rep(c("0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40"), nm)) 
    colnames(tmp)[6] <- "pval_casescontrol_m1" #to facilitate melt later
    
    tmp_melt <- tmp %>%
      pivot_longer(., cols = 2:3, names_to = c("group","remove"),
                   names_pattern = "([a-z]{5,7})(\\_rc\\_.+)",
                   values_to = paste0("rc_", i)) %>%
      mutate(., pval = ifelse(.$group == "cases", .$pval_cases_m1, .$pval_control_m1)) %>%
      select(-remove)
    assign(paste0("tmp_", i), tmp_melt)}
  
  dat_melt_comb <- cbind(dat_melt, tmp_ori, tmp_rand)
  dat_melt_comb <- dat_melt_comb[, !duplicated(colnames(dat_melt_comb))]
  
  final_melt <- dat_melt_comb
  
  #tidy up aesthetics 
  final_melt$group <- gsub("cases","Cases", final_melt$group)
  final_melt$group <- gsub("control","Controls", final_melt$group)
  final_melt$thresh <- as.factor(final_melt$thresh) %>% fct_relevel(.,"0.40", "0.35", "0.30", "0.25", "0.20", "0.15", "0.10")
  
  return(final_melt)
}

ukb_dat <- org_dat_FUN(grp_pval_cidi_fa)

plt_dat <- ukn_dat
plt_dat <- plt_dat %>% filter(thresh == "0.35") #best performing threshold

#f1: cases; f2: controls; f3: cases and controls
f1_dat <- plt_dat %>% filter(group == "Cases" & k < 63) %>%  
  mutate(sig_xloc = rowMeans(.[,c("xmin","xmax")], na.rm=T)) 
f2_dat <- plt_dat %>% filter(group == "Controls" & k < 63) %>% 
  mutate(sig_xloc = rowMeans(.[,c("xmin","xmax")], na.rm=T)) 
f3_dat <- plt_dat %>% filter(k < 63) 

#run plotly
#https://plotly.com/r/line-and-scatter/

f1 <- plot_ly(f1_dat, x = ~k, y = ~rc_ori, mode = "lines+markers",type = "scatter", 
              name = paste0("\u03A6","ori(k)"), color = I("#051006"), opacity = 0.5) %>%
  add_trace(y = ~rc_rand, mode = "lines+markers",  type = "scatter",
            name = paste0("\u03A6","rand(k)"),
            marker = list(color = 'grey'),
            line = list(color = 'grey'),
            opacity = 0.5) %>% 
  add_trace(x = ~k, y = ~rc_norm, mode = "lines+markers", type = "scatter", 
            marker = list(color = 'FF0000'), opacity = 1,
            line = list(color = 'FF0000'),
            yaxis = "y2", name = paste0("\u03A6","norm(k)")) %>%
  layout(yaxis2 = list(range = list(0.99, 1.20),
                       automargin = T, 
                       tickfont = list(color = "red"),
                       showline = T, linecolor = "red",
                       overlaying = "y", side = "right",
                       title = list(text = paste0("<b>\u03A6","norm(k)<b>"), 
                                    font = list(color = "red"))), #add second y axis
         title = '<b>Cases<b>',
         xaxis = list(title="<b>k<b>"), 
         yaxis = list(range = list(0.75, 1.05), title=paste0("<b>\u03A6","(k)<b>")), 
         #add rectangle to show significant region
         shapes = list(list(type = "rect",
                            fillcolor = "gray", line = list(color = "#c5d7dd"), 
                            opacity = 0.1, layer = "below",
                            x0 = f1_dat$xmin[1], x1 = f1_dat$xmax[1], xref = "x",
                            y0 = 0.75, y1 = f1_dat$yrect[1], yref = "y1"),
                       list(type = "line", x0 = 0, x1 = max(f1_dat$k),
                            y0 = 1, y1 = 1, yref = "y2",
                            line = list(color = "red", dash="dot"))),
         annotations = list(y = 1.025, x = f1_dat$sig_xloc[1], 
                            text = "<i>rich club organisation<i>", showarrow = F, 
                            yref = "y1"), 
         showlegend = F)

f2 <- plot_ly(f2_dat, x = ~k, y = ~rc_ori, mode = "lines+markers",type = "scatter", 
              name = paste0("\u03A6","ori(k)"), color = I("#051006"), opacity = 0.5) %>%
  add_trace(y = ~rc_rand, mode = "lines+markers",  type = "scatter",
            name = paste0("\u03A6","rand(k)"),
            marker = list(color = 'grey'),
            line = list(color = 'grey'),
            opacity = 0.5) %>% 
  add_trace(x = ~k, y = ~rc_norm, mode = "lines+markers", type = "scatter", 
            marker = list(color = '#0000CC'),
            line = list(color = '#0000CC'), opacity = 1,
            yaxis = "y2", name = paste0("\u03A6","norm(k)")) %>%
  layout(yaxis2 = list(range = list(0.99, 1.20),
                       automargin = T,
                       tickfont = list(color = "#0000CC"),
                       showline = T, linecolor = "#0000CC",
                       overlaying = "y", side = "right",
                       title = list(text = paste0("<b>\u03A6","norm(k)<b>"), 
                                    font = list(color = "#0000CC"))), #add second y axis
         title = '<b>Controls<b>',
         xaxis = list(title="<b>k<b>"), 
         yaxis = list(range = list(0.75, 1.05), title=paste0("<b>\u03A6","(k)<b>")), 
         #add rectangle to show significant region
         shapes = list(list(type = "rect",
                            fillcolor = "gray", line = list(color = "#c5d7dd"), 
                            opacity = 0.1, layer = "below",
                            x0 = f2_dat$xmin[1], x1 = f2_dat$xmax[1], xref = "x",
                            y0 = 0.75, y1 = f2_dat$yrect[1], yref = "y1"),
                       list(type = "line", x0 = 0, x1 = max(f1_dat$k),
                            y0 = 1, y1 = 1, yref = "y2",
                            line = list(color = "#0000CC", dash="dot"))),
         annotations = list(y = 1.025, x = f2_dat$sig_xloc[1], 
                            text = "<i>rich club organisation<i>", showarrow = F), 
         showlegend = F)


f3 <- plot_ly(f3_dat, x = ~k, y = ~rc_norm, mode = "lines+markers",type = "scatter",
              split = ~group, color = ~group, colors = c("#FF0000", '#0000CC'), opacity = 1) %>%
  layout(title = '<b>Cases VS Controls<b>',
         xaxis = list(title="<b>k<b>"), 
         yaxis = list(range = list(0.99, 1.20), title=paste0("<b>\u03A6","norm(k)<b>")),
         shapes = list(list(type = "line", x0 = 0, x1 = max(f3_dat$k),
                            y0 = 1, y1 = 1, yref = "y1",
                            line = list(color = "black", dash="dot"))))

#repeat for GS
