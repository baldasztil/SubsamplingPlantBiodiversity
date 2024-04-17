
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
library(vegan)
library(DEoptim)
library(doSNOW)

# Defining functions -----------------------------------------------------------

output <- function(species=NULL,best_corr=NULL){
  me <- data.frame(
    Species = species,
    best_corr = best_corr
  )
  
  return(me)
}

cor.gw <- function(rich_rel_shp) {
  
  bw <- try(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
                   adaptive = T, 
                   longlat = T))
  
  if (!is.numeric(bw)) {
    bw <- bw_global
  }
  
  stats <- try(as.data.frame(gwss(rich_rel_shp, rich_rel_shp, 
                                  vars = c("richness", "richness_sample"),
                                  adaptive = T, bw = bw, longlat = T)$SDF))
  
  if (class(stats)[1]=="try-error") {
    error_data <- 10000
    return(error_data)
  }
  
  else {
    stats_total <-  cbind(stats, rich_rel_shp) %>% 
      dplyr::select(cor.sp_gwr = Spearman_rho_richness.richness_sample,
                    LEVEL3_COD, 
                    continent = LEVEL1_NAM) %>% 
      replace(is.na(.), 0) %>% 
      setDT()
    
    funs <- lst(mean, sd)
    
    global <- stats_total %>%  
      summarise(across(
        .cols = where(is.numeric), 
        .fns = list(mean = mean, sd = sd, med = median),  
        .names = "{col}_{fn}"
      )) %>% 
      mutate(
             dim = unique(rich_rel_shp$dim)) 
    return(global)
  }
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
  prune_tree <- try(match_phylo_comm(output_tree,subset_matrix))
  
  if (class(prune_tree)[1] =="try-error")  {
    global_er <- Inf 
    return(global_er)
  } else 
    
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
  
  rich_rel_shp <- lapply(split(rich_rel,rich_rel$dim), FUN = function(x)
    as(x, "Spatial"))
  
  
  correlation <- try(lapply(rich_rel_shp, cor.gw) %>% 
    rbindlist())
  
  if (class(correlation)[1] =="try-error")  {
    global_er <- Inf 
    return(global_er)
  } else 
  
  cor.gwr <- sum(correlation$cor.sp_gwr_mean, na.rm = T)
  
  if (cor.gwr > max(maxima)) {
    assign("maxima", cor.gwr, envir = .GlobalEnv)
    #maxima <- list(maxima, as.list(cor.gwr))
    ids <- unique(species_sample$plant_name_id)
    write.table(c(cor.gwr,ids), paste0("./data/",
                                       length(ids),"_target_population_corgw_gbif.txt"))
  } 
  
  return(-cor.gwr)
}

deoptim_lapply <- function(sp_num){
  assign("maxima", 0, envir = .GlobalEnv)
  de_result <- DEoptim(subsampling.plants.solo, lower = sp_num, upper = sp_num,
                       DEoptim.control(VTR = -Inf, strategy = 1, NP = 10, 
                                       bs = FALSE, itermax = 100, CR = 0.9, F = 0.8, 
                                       storepopfreq = 1, trace = T, c = 0.55, 
                                       steptol = 100))
  return(output(de_result$optim$bestmem, de_result$optim$bestval*-1))
}


# Import data ------------------------------------------------------------------
gbif_dist <- fread("data/dist_gbif.txt")  %>% 
  rename(area_code_l3 = LEVEL3_COD)
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

gbif_dist_calc <- gbif_dist %>% 
  left_join(plants_full, by = "plant_name_id")

table_data <- table(gbif_dist_calc$area_code_l3, gbif_dist_calc$growth_form)


shannon_div <- as.data.frame(diversity(table_data, index = "shannon")) 
shannon_div$C <- rownames(shannon_div)
names(shannon_div) <- c("richness_fun", "LEVEL3_COD")

richness_patterns_fun <- shannon_div 

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

plantlist_dist_phylo_growth <- gbif_dist %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form, family) %>% 
  mutate(label = gsub(" ", "_", taxon_name))

dist.mat <- long2sparse(plantlist_dist_phylo_growth, grids = "LEVEL3_COD", species = "label")
prune_tree <- match_phylo_comm(output_tree,dist.mat)


richness_patterns_phy <- as.data.frame(PD(dist.mat, prune_tree$phy), 
                                       col.names =c("richness")) %>% 
  mutate(LEVEL3_COD = rownames(.)) %>% 
  rename(richness_phy = `PD(dist.mat, prune_tree$phy)`) 


richness_patterns_rich <- gbif_dist_calc %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = length(unique(plant_name_id))) %>% 
  rename(LEVEL3_COD = area_code_l3)

richness_patterns_tax <- gbif_dist_calc %>% 
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

output_model <- lapply(rich_overall_bru_mid_test, cor.gw) %>% 
  rbindlist() %>% 
  summarise(sum(cor.sp_gwr_mean))




# compare models ---------------------------------------------------------------

start <- Sys.time()
deoptim_lapply(50000)
Sys.time() - start

set.seed(1234)

solutions_listed <- mclapply(seq(2, nrow(plantlist_names) / 100, 1), mc.cores = 70, 
                             deoptim_lapply)

solutions <- do.call(rbind, solutions_listed)
write.table(solutions, "./data/optimised/DEoptim_results_corgw_gbif.txt")
