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
  `Wu-1` = "#393b79",
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
  "../../data/01_sars_cov_2_antigen_colors.RDS"
)

# make coords ------------------------------------------------------------------

example_map = Racmacs::read.acmap("../data/initial_example_map.ace")

initial_ag_coords = agCoords(example_map)[,1:2]

initial_ag_coords["XBB.1.5",] = initial_ag_coords["XBB.1.5",] + c(0.8, 0.2)

initial_ag_coords[c("JN.1", "KP.2.3", "KP.3.1.1"),] = initial_ag_coords[
  c("JN.1", "KP.2.3", "KP.3.1.1"),] + matrix(c(0.8, 0.2), 3, 2, T)

full_ag_coords = rbind(
  initial_ag_coords,
  
  B.1.621 = initial_ag_coords["Wu-1",] + c(1, 1.5),
  
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
rownames(full_sr_coords) = paste0(rownames(full_sr_coords), "|", 1:3)
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

log_peak_titers = log2(2560) + rnorm(ncol(distances), 0, 0.8)

log_titers = array(
  rep(log_peak_titers, each = nrow(distances)),
  dim = dim(distances)
) - distances

titers = 2**log_titers

write.csv(
  titers,
  file = "../../data/01_sars_cov_2_titerdata.csv"
)

# sequence data ----------------------------------------------------------------

GISAID_ids = c(
  `Wu-1` = "EPI_ISL_402124",
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
aligned_spike_reference_dna = seqinr::read.fasta("../data/wu1_reference.fasta")[[2]]

# aligning
full_genomes = seqinr::read.fasta("../data/gisaid_hcov-19_2025_06_17_12.fasta")

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
  "../../data/01_spike_sequences.fasta"
  )

# test -------------------------------------------------------------------------
rm(list = ls())

titers = read.csv("../../data/01_sars_cov_2_titerdata.csv", row.names = 1)

reconstructed_map = Racmacs::acmap(titer_table = titers) %>%
  Racmacs::optimizeMap(
    number_of_dimensions = 2,
    number_of_optimizations = 1000
  )

antigen_colors = readRDS("../../data/01_sars_cov_2_antigen_colors.RDS")

Racmacs::agFill(reconstructed_map) = antigen_colors[agNames(reconstructed_map)]
Racmacs::srOpacity(reconstructed_map) = 0.4

