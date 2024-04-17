
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(parallel)
library(doParallel)
library(sf)
library(GWmodel)    # to undertake the GWR
library(foreach)
library(future.apply)
library(rWCVP)



# Defining functions -----------------------------------------------------------
std.error <- function(x) sd(x)/sqrt(length(x))
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  
  plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist)
  if (nrow(plantlist_names_left) > spec_n) {
    species_sample <- sample_n(plantlist_names_left, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names_left
  }
  cumulative_namelist <- c(cumulative_namelist, species_sample$plant_name_id)
  species <-plantlist_names %>% 
    filter(plant_name_id %in%  cumulative_namelist)
  
  dist <- dist_native %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # richness patterns across brus
  rich_rel <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL3_COD = area_code_l3) %>% 
    right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist),
           richness_sample = sqrt(richness_sample),
           richness = sqrt(richness)) %>% 
    st_sf() 
    
    rich_rel_shp <- rich_rel %>% 
    as("Spatial")
  

  m_spatial <- spmodel::splm(richness ~ richness_sample, data=rich_rel,
                      spcov_type =  c("exponential", "spherical", "gaussian", 
                                      "triangular", "circular", "cubic", 
                                      "pentaspherical", "cosine", "wave", 
                                      "jbessel", "gravity", "rquad",
                                      "magnetic", "matern", "cauchy", "pexponential", 
                                      "none"))
  model_comp <- glances(m_spatial)
  model_comp$sp <- length(species$plant_name_id)
  return(model_comp)
}



# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns <- fread("data/richness_patterns.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")


# manipulate data --------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id,continent_code_l1,area_code_l3)

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)


rich_overall_bru <- richness_patterns

rich_overall_bru_mid <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") 

rich_overall_bru_shp <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") 


# compare models ---------------------------------------------------------------

plan(multisession, workers = 6)
modelcomp_test <- future_lapply(sample(1:length(plantlist_names$plant_name_id), 100), 
                                subsampling.plants, future.seed = T)

output <- rbindlist(modelcomp_test)

out_stats <- output %>% 
  group_by(sp) %>% 
  summarise_all(funs(n(),mean,median,sd, std.error))

out_first  <- output %>% 
  slice_min(n= 1, sp)


output_temp <- output %>% 
  dplyr::select(model) %>% 
  summarise(model = unique(model))

best_performer_loglike <- output %>% 
  group_by(sp) %>% 
  slice_max(logLik, n = 3) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  tally(name = "loglike")

best_performer_AICc <- output %>% 
  group_by(sp) %>% 
  slice_min(AICc, n = 3) %>% 
  group_by(model) %>% 
  tally(name = "AICc")

best_performer_pseudoR <- output %>% 
  group_by(sp) %>% 
  slice_max(pseudo.r.squared, n = 3) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  tally(name = "pseudoR")

best_performer_dev <- output %>% 
  group_by(sp) %>% 
  slice_min(deviance, n = 3) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  tally(name = "deviance")

output_comb <- output_temp %>% 
  left_join(best_performer_loglike, by = "model") %>% 
  left_join(best_performer_AICc, by = "model") %>% 
  left_join(best_performer_pseudoR, by = "model") %>% 
  left_join(best_performer_dev, by = "model") %>% 
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = T))

  
