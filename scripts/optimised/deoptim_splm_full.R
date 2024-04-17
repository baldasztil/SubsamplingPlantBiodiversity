
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(spmodel)
library(ape)
library(phyloregion)
library(GWmodel)
library(doSNOW)
library(vegan)
library(DEoptim)

# Defining functions -----------------------------------------------------------

output <- function(species=NULL,best_corr=NULL){
  me <- data.frame(
    Species = species,
    best_corr = best_corr
  )
}

spatial_lm <- function(x)  {
  model <- splm(richness ~ richness_sample, data = x,  
                spcov_type =  c("exponential"))
  model$dim <- x$dim
  return(model)
}

model_preds <-  function(model) {
  temp <- as.data.table(model$obdata)[, c("LEVEL1_NAM", "LEVEL3_COD", "richness", "richness_sample")] %>%
    mutate(
           residuals = residuals(model, "pearson"), 
           dim = unique(model$dim)) 
  return(temp)
}

glance_splm <-  function(model) {
  temp <-  glance(model) %>% 
    mutate(continent = "GLOBAL", 
           model_id = 1, 
           coeff = model$coefficients$fixed[2],
           coeffse = sqrt(diag(vcov(model)))[2], 
           dim = unique(model$dim))
  
  return(temp)
}

subsampling.plants.solo <- function(spec_n) {
  
  if (nrow(plantlist_names) > spec_n) {
    species_sample <- sample_n(plantlist_names, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names
  }
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  # phylo richness
  subset_matrix <- dist.mat[ ,colnames(dist.mat) %in% dist$label]
  prune_tree <- match_phylo_comm(output_tree,subset_matrix)
  
  # doesnt work right now with old tree
  rich_rel_pd <- as.data.frame(PD(subset_matrix, prune_tree$phy)) %>% 
    mutate(LEVEL3_COD = rownames(.)) %>% 
    rename(richness_phy = `PD(subset_matrix, prune_tree$phy)`)
  
  table_data_sample <- table(dist$LEVEL3_COD, dist$growth_form)
  rich_rel_fun <- as.data.frame(diversity(table_data_sample, 
                                          index = "shannon")) 
  rich_rel_fun$C <- rownames(rich_rel_fun)
  names(rich_rel_fun) <- c("richness_fun", "LEVEL3_COD")
  
  rich_rel_rich_tax <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness = length(unique(plant_name_id)), 
              richness_tax = length(unique(family)))
  
  richness_rel <- rich_rel_rich_tax %>% 
    left_join(rich_rel_pd, by = "LEVEL3_COD") %>% 
    left_join(rich_rel_fun, by = "LEVEL3_COD")
  
  
  rich_rel_long <- richness_rel %>% 
    pivot_longer(cols = starts_with("richness"),
                 names_to = "dim",
                 values_to = "richness_sample") %>% 
    mutate(richness_sample_scaled = scale(richness_sample)[,1]) 
  
  rich_rel <-  rich_rel_long %>% 
    right_join(rich_overall_bru, by = c("LEVEL3_COD", "dim")) %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample),
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf() 
  
  
  
  model <- try(lapply(split(rich_rel,rich_rel$dim),spatial_lm))
  
  if (class(model)[1] =="try-error")  {
    global_er <- 10000
    return(global_er)
  } else 
    
    global_mod <- lapply(model, glance_splm) %>% 
    rbindlist() %>% 
    mutate(sp = length(species_sample$plant_name_id))
  
  
  #rich_rel_preds  <- lapply(model, model_preds) %>% 
   # rbindlist() %>% 
    #group_by(dim) %>% 
    #summarise(resids = sum(residuals)^2)
  
  #global_mod$resids <- rich_rel_preds$resids
  splm.gwr <- sum(global_mod$pseudo.r.squared, na.rm = T)
  
  if (splm.gwr > max(maxima)) {
    assign("maxima", splm.gwr, envir = .GlobalEnv)
    #maxima <- list(maxima, as.list(splm.gwr))
    ids <- unique(species_sample$plant_name_id)
    write.table(c(splm.gwr,ids), paste0(
      "./data/optimised/pops/",
      length(ids),"_target_population_splm.txt"))
  } 
  return(-splm.gwr)
}

deoptim_lapply <- function(sp_num){
  assign("maxima", 0, envir = .GlobalEnv)
  de_result <- DEoptim(subsampling.plants.solo, lower = sp_num, upper = sp_num,
                       DEoptim.control(VTR = -Inf, strategy = 2,
                                       bs = FALSE, itermax = 10, CR = 0.5, F = 0.8, 
                                       storepopfreq = 1, trace = T, c = 0.55, 
                                       steptol = 100))
  return(output(de_result$optim$bestmem, de_result$optim$bestval*-1))
}


# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full_raw <- fread("data/wcvp_accepted_merged.txt")

output_tree <- read.tree("data/output_tree_20012024.tre")

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

growth <- read.csv("data/growths_rf.csv") %>% 
  dplyr::select(growth_form = apd_plant_growth_form, plant_name_id)

plants_full <- plants_full_raw %>% 
  left_join(growth, by = "plant_name_id")

dist_native_calc <- dist_native %>% 
  left_join(plants_full, by = "plant_name_id")

table_data <- table(dist_native_calc$area_code_l3, dist_native_calc$growth_form)


shannon_div <- as.data.frame(diversity(table_data, index = "shannon")) 
shannon_div$C <- rownames(shannon_div)
names(shannon_div) <- c("richness_fun", "LEVEL3_COD")

richness_patterns_fun <- shannon_div 

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

plantlist_dist_phylo_growth <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form,
                family) %>% 
  mutate(label = gsub(" ", "_", taxon_name))

dist.mat <- long2sparse(plantlist_dist_phylo_growth, grids = "LEVEL3_COD", species = "label")
prune_tree <- match_phylo_comm(output_tree,dist.mat)


richness_patterns_phy <- as.data.frame(PD(dist.mat, prune_tree$phy), 
                                       col.names =c("richness")) %>% 
  mutate(LEVEL3_COD = rownames(.)) %>% 
  rename(richness_phy = `PD(dist.mat, prune_tree$phy)`) 


richness_patterns_rich <- dist_native_calc %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = length(unique(plant_name_id))) %>% 
  rename(LEVEL3_COD = area_code_l3)

richness_patterns_tax <- dist_native_calc %>% 
  group_by(area_code_l3) %>% 
  summarise(richness_tax = length(unique(family))) %>% 
  rename(LEVEL3_COD = area_code_l3)


richness_patterns <-richness_patterns_rich %>% 
  left_join(richness_patterns_tax, by = "LEVEL3_COD") %>% 
  left_join(richness_patterns_phy, by = "LEVEL3_COD") %>% 
  left_join(richness_patterns_fun, by = "LEVEL3_COD") %>% 
  left_join(continent_names, by = "LEVEL3_COD") 

rich_patterns_long <- richness_patterns %>% 
  pivot_longer(cols = starts_with("richness"),
               names_to = "dim",
               values_to = "richness") %>% 
  mutate(richness_scaled = scale(richness)[,1])


rich_overall_bru <- rich_patterns_long


rich_overall_bru_mid <- rich_patterns_long %>% 
  left_join(midpoints_red, by = "LEVEL3_COD")  %>% 
  mutate(richness2 = richness) %>% 
  st_as_sf()

rich_overall_bru_shp <- rich_patterns_long %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_overall_bru_mid_shp <- rich_overall_bru_mid %>% 
  mutate(richness_sample = richness) %>% 
  as("Spatial")

rich_overall_bru_mid_test <- split(rich_overall_bru_mid_shp,rich_overall_bru_mid_shp$dim)


bw_global <-bw.gwr(richness ~ richness2, rich_overall_bru_mid_shp, 
                   approach = "AIC",
                   adaptive = T)


setDT(plantlist_names)
setDT(plantlist_dist_phylo_growth)
setDT(rich_overall_bru)
setDT(midpoints_red)


deoptim_lapply(2)

# compare models ---------------------------------------------------------------

library(profvis)
profvis(deoptim_lapply(2))
solutions_listed <- mclapply(seq(2, ceiling(nrow(plantlist_names) / 100) , 1), mc.cores = 10, deoptim_lapply)

solutions <- do.call(rbind, solutions_listed)
write.table(solutions, "./data/optimised/DEoptim_results_splm.txt")


