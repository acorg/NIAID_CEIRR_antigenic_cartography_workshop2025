rm(list = ls())

library(Racmacs)
library(tidyverse)

set.seed(1)

# antigen colors ---------------------------------------------------------------
og_map <- read.acmap("data/maps/01_sars_cov_2_map.ace")

map <- subsetMap(og_map, antigens = c("Wu1", "B.1.621", "BA.2", "BA.5", "XBB.1.5", "HK.3", "JN.1", "KP.2.3"))
map <- realignMap(optimizeMap(map, 2, 1000), og_map)

srGroups(map) <- sapply(srNames(map), function(x){
  temp_split <- strsplit(x, "\\_")[[1]]
  paste(temp_split[1:(length(temp_split)-1)], collapse = ".")
  
})


map <- subsetMap(map, sera = srNames(map)[srGroups(map) %in% agNames(map)])
map <- realignMap(optimizeMap(map, 2, 1000), og_map)

srOutline(map) <- agFill(map)[as.character(srGroups(map))]

# base map for which we will add multi exposure sera (condensation)
plot(map)

# create artificial multi exposure sera
reduce_serum_fc <- function(single_map, sr_target = "Wu1", target_slope = 0.8, sr_new = "BT", logtiter_increase = 0){
  
  single_wu1 <- logtiterTable(single_map)[,grepl(sr_target, colnames(titerTable(single_map)))] + logtiter_increase
  single_table_dist <- tableDistances(single_map)[,grepl(sr_target, colnames(titerTable(single_map)))]
  
  single_table_dist <- apply(single_table_dist, 1:2, function(x){
    as.numeric(gsub(">", "", x))
  })
  
  delta_bt <- abs(sweep(target_slope*apply(single_table_dist, c(1,2), as.numeric), MARGIN = 2, single_wu1[1,c(1:3)], `-`))
  dimnames(delta_bt) <- dimnames(single_wu1)
  colnames(delta_bt) <- gsub(sr_target, sr_new, colnames(delta_bt))
  delta_bt <- apply(delta_bt, 1:2, function(x)2^x*10)

  return(delta_bt)
}

# add higher titers for reexposed Ags
wt_xbb15 <- reduce_serum_fc(map, "Wu1", target_slope = 0.3, sr_new = "Unknown 1", logtiter_increase = 1.5)
wt_xbb15[c("XBB.1.5", "HK.3"), ] <- wt_xbb15[c("XBB.1.5", "HK.3"), ]*1.2
xbb15_reinf <- reduce_serum_fc(map, "BA.2", target_slope = 0.4, sr_new = "Unknown 2", logtiter_increase = 2)
xbb15_reinf[c("XBB.1.5", "HK.3", "JN.1", "KP.2.3"), ] <- xbb15_reinf[c("XBB.1.5", "HK.3", "JN.1", "KP.2.3"), ]*1.4
sample3 <- reduce_serum_fc(map, "B.1.621", target_slope = 0.2, sr_new = "Unknown 3", logtiter_increase = 1.5)
sample3[c("XBB.1.5", "HK.3", "JN.1"), ] <- sample3[c("XBB.1.5", "HK.3", "JN.1"), ]*1.2

# make new map
new_tt <- cbind(titerTable(map), wt_xbb15, xbb15_reinf, sample3)
new_map <- make.acmap(new_tt)
srGroups(new_map) <- sapply(srNames(new_map), function(x){
  temp_split <- strsplit(x, "\\_")[[1]]
  paste(temp_split[1:(length(temp_split)-1)], collapse = ".")
})

plot(procrustesMap(new_map, map))
## check new map
# agFill(new_map) <- agFill(map)
# srOutline(new_map) <- c(srOutline(map), rep("grey", 3), rep("grey30", 3))

save.acmap(new_map, "data/maps/03_troubleshooting_map2a.ace")

# ------------------------ make map below of very early sampled sera, e.g. close 
# -------------------------to detection threshold for all except 1 antigen.
create_early_sampled_sera <- function(single_map, sr_target = "Wu1", target_lod = log2(40/10), sr_new = "BT", target_sd = 2,
                                      homologous_ag = "Wu1"){
  
  single_wu1 <- logtiterTable(single_map)[,grepl(sr_target, colnames(titerTable(single_map)))]
  delta_bt <-matrix(target_lod + rnorm(length(single_wu1), mean = 0, sd = target_sd),
                    nrow = nrow(single_wu1),
                    ncol = ncol(single_wu1))
  dimnames(delta_bt) <- dimnames(single_wu1)
  delta_bt[homologous_ag, ] <- target_lod + 1 + rnorm(ncol(single_wu1), mean = 0, sd = 0.7)
  
  colnames(delta_bt) <- gsub(sr_target, sr_new, colnames(delta_bt))
  delta_bt <- apply(delta_bt, 1:2, function(x){
    temp <- 2^x*10
    if(temp < 2^target_lod*10){
      paste0("<", 2^target_lod*10)
    } else {
      paste0(temp)
    }
    })
  
  return(delta_bt)
}

# simulate early sera with random noise around detection limit, but higher titers against their homolgous ag
early_wu1 <- create_early_sampled_sera(map, "Wu1", log2(40/10), sr_new = "Wu1 Day 4", target_sd = 1, homologous_ag = "Wu1")
early_jn1 <- create_early_sampled_sera(map, "JN.1", log2(40/10), sr_new = "JN.1 Day 4", target_sd = 1, homologous_ag = "JN.1")
early_ba2 <- create_early_sampled_sera(map, "BA.2", log2(40/10), sr_new = "BA.2 Day 4", target_sd = 1, homologous_ag = "BA.2")

# make new map
new_tt <- cbind(titerTable(map), early_wu1, early_jn1, early_ba2)

new_tt <- apply(new_tt, 1:2, function(x){
  temp <- as.numeric(x)
  if(!is.na(temp)){
    if(temp < 40){
      "<40"
    } else {
      x
    }
  } else {
    x
  }
})


new_map <- make.acmap(new_tt)
srGroups(new_map) <- sapply(srNames(new_map), function(x){
  temp_split <- strsplit(x, "\\_")[[1]]
  paste(temp_split[1:(length(temp_split)-1)], collapse = ".")
})

save.acmap(new_map, "data/maps/03_troubleshooting_map2b.ace")

# #check new map
# srOutline(new_map) <- c(srOutline(map), rep("grey", 3), rep("grey30", 3), rep("grey80", 3))
# agFill(new_map) <- agFill(map)
# plot(new_map, plot_stress = TRUE)
# plot(procrustesMap(map, new_map))
# RacViewer(new_map)

