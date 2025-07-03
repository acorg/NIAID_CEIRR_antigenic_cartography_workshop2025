setwd(here::here("internal/code/"))

rm(list = ls())

library(Racmacs)
library(magrittr)
library(stringr)
library(purrr)
library(seqUtils) # devtools::install_github("samt123/seqUtils")

set.seed(1)

# antigen colors ---------------------------------------------------------------

antigen_colors = c(
  `Wu1` = "#393b79",
  B.1.621 = "#e7ba52",
  
  BA.2 = "#5B004C",
  BA.2.75 = "#95004b",
  BA.2.12.1 = "#8F5180",
  
  BA.5 = "#e31d5e",
  BQ.1.1 = "#ec6391",
  BF.7 ="#FF6A8B",
  
  XBB.1.5 = "#ff6500",
  HV.1 = "#FF9A81",
  HK.3 = "#FF815A",
  
  JN.1 = "#7ac900",
  KP.2.3 = "#BDE563",
  KP.3.1.1 = "#C3CF1F"
)

saveRDS(
  antigen_colors,
  "../../data/other/01_sars_cov_2_antigen_colors.RDS"
)

# make coords ------------------------------------------------------------------

example_map = Racmacs::read.acmap("../data/initial_example_map.ace")
agNames(example_map)[agNames(example_map) == "Wu-1"] = "Wu1"

initial_ag_coords = agCoords(example_map)[,1:2]

initial_ag_coords["XBB.1.5",] = initial_ag_coords["XBB.1.5",] + c(0.8, 0.2)

initial_ag_coords[c("JN.1", "KP.2.3", "KP.3.1.1"),] = initial_ag_coords[
  c("JN.1", "KP.2.3", "KP.3.1.1"),] + matrix(c(0.8, 0.2), 3, 2, T)

full_ag_coords = rbind(
  initial_ag_coords,
  
  B.1.621 = initial_ag_coords["Wu1",] + c(1, 1.5),
  
  BQ.1.1 = initial_ag_coords["BA.5",] + c(0.3, 0.35),
  BF.7 = initial_ag_coords["BA.5",] + c(0.15, -0.3),
  
  BA.2.75 = initial_ag_coords["BA.2",] + c(0.5, -0.4),
  BA.2.12.1 = initial_ag_coords["BA.2",] + c(0.1, 0.2),
  
  HV.1 = initial_ag_coords["XBB.1.5",] + c(0.2, -0.3),
  HK.3 = initial_ag_coords["XBB.1.5",] + c(0.3, -0.4)
)

full_ag_coords = full_ag_coords[names(antigen_colors),]

plot(full_ag_coords)
text(full_ag_coords, labels = rownames(full_ag_coords))


full_sr_coords = full_ag_coords[rep(seq_len(nrow(full_ag_coords)), each = 3),]
rownames(full_sr_coords) = paste0(rownames(full_sr_coords), "_", 1:3)
full_sr_coords = full_sr_coords + rnorm(length(full_sr_coords), 0, 0.5)

# construct titers -------------------------------------------------------------

# empty distance table
distances = matrix(
  NA,
  nrow = nrow(full_ag_coords),
  ncol = nrow(full_sr_coords),
  dimnames = list(
    rownames(full_ag_coords),
    rownames(full_sr_coords)
  )
)

for (i in seq_len(nrow(full_ag_coords))){
  for (j in seq_len(nrow(full_sr_coords))){
    distances[i, j] = sqrt(sum((full_ag_coords[i, ] - full_sr_coords[j,])**2))
  }
}

distances_serumnoise = distances * 
  matrix(2**rnorm(ncol(distances), 0, 0.2), nrow(distances), ncol(distances), byrow = T)

log_peak_titers = log2(1280) + rnorm(ncol(distances_serumnoise), 0, 1)

log_titers = array(
  rep(log_peak_titers, each = nrow(distances_serumnoise)),
  dim = dim(distances_serumnoise)
) - distances_serumnoise

log_titers = log_titers + rnorm(length(log_titers), 0, 0.7)

titers = 2**log_titers
titers[titers < 40] = "<40"

# subsetting and saving
serum_variants = colnames(titers) %>%
  str_split(fixed("_")) %>%
  map_chr(1)

## square map 
write.csv(
  titers,
  file = here::here("data", "titerdata", "01_sars_cov_2_titerdata_square.csv")
)


## standard map
standard_map_serum_variants = c("Wu1", "B.1.621", "BA.2", "BA.2.75", "BA.5", "BQ.1.1", "XBB.1.5", "JN.1", "KP.3.1.1")
titers_standard = titers[, serum_variants %in% standard_map_serum_variants]
write.csv(
  titers_standard,
  file = here::here("data", "titerdata", "01_sars_cov_2_titerdata.csv")
)

## missing map
missing_map_serum_variants = standard_map_serum_variants[
  !standard_map_serum_variants %in% c("XBB.1.5", "JN.1", "KP.3.1.1")]

titers_missing = titers[, serum_variants %in% missing_map_serum_variants]
write.csv(
  titers_missing,
  file = here::here("data", "titerdata", "02_sars_cov_2_titerdata_missing.csv")
)

## with two extra early sera
extraearly_map_serum_variants = c(missing_map_serum_variants, "BF.7", "BA.2.12.1")
titers_extraearly = titers[, serum_variants %in% extraearly_map_serum_variants]
write.csv(
  titers_extraearly,
  file = here::here("data", "titerdata", "02_sars_cov_2_titerdata_missing_plus_BA275_BQ11.csv")
)


## with extra JN.1 and XBB.1.5 sera
extralate_map_serum_variants = c(missing_map_serum_variants, "XBB.1.5", "JN.1")
titers_extralate = titers[, serum_variants %in% extralate_map_serum_variants]
write.csv(
  titers_extralate,
  file = here::here("data", "titerdata", "02_sars_cov_2_titerdata_missing_plus_XBB15_JN1.csv")
)

## no triplicate

titers_not_trip = titers[, str_ends(colnames(titers), fixed("_1"))]
write.csv(
  titers_not_trip,
  file = here::here("data", "titerdata", "02_sars_cov_2_titerdata_no_triplicate.csv")
)

# discordant map
# distances_discordant = distances
# distances_discordant[, serum_variants %in% c("JN.1")] = 
#   distances_discordant[, serum_variants %in% c("JN.1")] * 3
# 
# log_peak_titers_discordant = log2(2560) + rnorm(ncol(distances_discordant), 0, 1)
# 
# log_titers_discordant = array(
#   rep(log_peak_titers_discordant, each = nrow(distances_discordant)),
#   dim = dim(distances_discordant)
# ) - distances_discordant
# 
# log_titers_discordant = log_titers_discordant + rnorm(length(log_titers_discordant), 0, 0.7)
# 
# titers_discordant = 2**log_titers_discordant
# titers_discordant[titers_discordant < 40] = "<40"
# 
# titers_discordant = titers_discordant[, serum_variants %in% c("Wu1", "BA.2", "XBB.1.5", "JN.1")]
# 
# write.csv(
#   titers_discordant,
#   file = "../../data/01_sars_cov_2_titerdata_discordant.csv"
# )


# sequence data ----------------------------------------------------------------

GISAID_ids = c(
  `Wu1` = "EPI_ISL_402124",
  B.1.617.2 = "EPI_ISL_2657324",
  B.1.621 = "EPI_ISL_2828019",
  
  BA.2 = "EPI_ISL_7190366",
  BA.2.75 = "EPI_ISL_13583301",
  BA.2.12.1 = "EPI_ISL_10783322",
  
  BA.5 = "EPI_ISL_11207535",
  BQ.1.1 = "EPI_ISL_14752457",
  BF.7 = "EPI_ISL_15157485",
  
  XBB.1.5 = "EPI_ISL_16343798",
  HV.1 = "EPI_ISL_18271598",
  HK.3 = "EPI_ISL_17987533",
  
  JN.1 = "EPI_ISL_18300149",
  KP.2.3 = "EPI_ISL_18999330",
  KP.3.1.1 = "EPI_ISL_19176351"
)

cat(GISAID_ids, sep = "\n") # paste into GISAID search

# reference sequence
aligned_spike_reference_dna = seqinr::read.fasta(
  here::here("internal", "data", "wu1_reference.fasta")
  )[[2]]

# aligning
full_genomes = seqinr::read.fasta(
  here::here("internal", "data", "gisaid_hcov-19_2025_06_17_12.fasta")
)

aligned_spike_dnas = mafft_align(
  full_genomes %>% map_chr(paste, collapse = "") %>% map_chr(toupper),
  paste(aligned_spike_reference_dna, collapse = "")
)

names(aligned_spike_dnas) = setNames(
  names(GISAID_ids),
  unname(GISAID_ids))[
    names(aligned_spike_dnas) %>%
      str_split(fixed("|")) %>%
      map_chr(2)]

aligned_spike_dnas = aligned_spike_dnas[names(antigen_colors)]

# write to file
seqinr::write.fasta(
  stringr::str_split(aligned_spike_dnas, ""),
  names(aligned_spike_dnas),
  here::here("data", "other", "01_spike_sequences.fasta")
  )

# test -------------------------------------------------------------------------
rm(list = ls())
source(here::here("code", "helpful_functions", "rotateMapCorrectly.R"))

saveMap = function(titer_path, optimize = T){
  
  map_path = str_replace(titer_path, "titerdata", "map") %>%
    fs::path_ext_remove()
  
  if (!optimize){
    map_path = paste0(map_path, "_not_opt")
  }
  
  map_path = fs::path_ext_set(map_path, ".ace")
  
  titers = read.csv(
    here::here("data", "titerdata", titer_path),
    row.names = 1,
  
  )
  map = Racmacs::acmap(titer_table = titers)
  
  antigen_colors = readRDS(here::here("data", "other", "01_sars_cov_2_antigen_colors.RDS"))
  
  Racmacs::agFill(map) = antigen_colors[agNames(map)]
  Racmacs::srOpacity(map) = 0.4
  
  if (optimize){
    map = Racmacs::optimizeMap(
      map,
      number_of_dimensions = 2,
      number_of_optimizations = 1000
    )
    
    map = rotateMapCorrectly(map)
  }
  
  Racmacs::save.acmap(
    map,
    here::here("data", "maps", map_path)
  )  
  
}

saveMap(
  "01_sars_cov_2_titerdata_square.csv",
  T
)


saveMap(
  "01_sars_cov_2_titerdata.csv",
  T
)

saveMap(
  "01_sars_cov_2_titerdata.csv",
  F
)

saveMap(
  "02_sars_cov_2_titerdata_missing.csv",
  T
)

saveMap(
  "02_sars_cov_2_titerdata_missing_plus_BA275_BQ11.csv",
  T
)

saveMap(
  "02_sars_cov_2_titerdata_missing_plus_XBB15_JN1.csv",
  T
)

saveMap(
  "02_sars_cov_2_titerdata_no_triplicate.csv",
  T
)