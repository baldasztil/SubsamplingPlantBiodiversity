
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
  if (nrow(plantlist_names) > spec_n) {
    species_sample <- sample_n(plantlist_names, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names
  }
  
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  # richness patterns across brus
  rich_rel_shp <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample)) %>% 
    st_sf() %>% 
    as("Spatial")
  
 
 bw <- try(bw.gwr(richness ~ richness_sample, data=rich_rel_shp,
               adaptive = T, 
               longlat = T))
  
  if (!is.numeric(bw)) {
    bw <- bw_global
  }
 
 stats <- try(as.data.frame(gwss(rich_rel_shp, rich_rel_shp, 
                                 vars = c("richness", "richness_sample"),
                                 adaptive = T, bw = bw, longlat = T)$SDF)) 

    if (class(stats)=="try-error") {
    error_data <- output_model %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate(id = output_model$id, 
             sp = nrow(species_sample))
    return(error_data)
  }
  
  else {
    
    stats_total <-  cbind(stats, rich_rel_shp) %>% 
      dplyr::select(cor.sp_gwr =Spearman_rho_richness.richness_sample,
                    cor.pe_gwr = Corr_richness.richness_sample, LEVEL3_COD, LEVEL1_NAM) %>% 
      replace(is.na(.), 0)
    
    funs <- lst(mean, sd)
    
    global <- stats_total %>%  
      mutate(id = LEVEL1_NAM) %>%  
      summarise(across(
        .cols = where(is.numeric), 
        .fns = list(mean = mean, sd = sd,  med = median),  
        .names = "{col}_{fn}"
      )) %>% 
      mutate(id = "GLOBAL") 
      
   continents <- stats_total %>% 
     group_by(LEVEL1_NAM) %>% 
     rename(id = LEVEL1_NAM) %>% 
     summarise(across(
       .cols = where(is.numeric), 
       .fns = list(mean = mean, sd = sd, med = median), 
       .names = "{col}_{fn}"
     )) 
   
   output_file <- rbind(continents, global) %>% 
     mutate(sp = nrow(species_sample)) %>% 
     mutate_if(is.numeric,
               round,
               digits = 4)

  # cumulative pattern
    if (nrow(species_sample) %in% seq(1,nrow(plantlist_names),20000)){
    a <-  paste0("There are ", nrow(species_sample) ," species in the subsample")
    write.table(a, paste0("data/fullsamples_test/checkpoints/Check@_future_sampling",
                          nrow(species_sample),"sp.txt")) 
  }
  
  return(output_file)
  }
}
write.sample <- function(x) {
  xx <- as.data.frame(samples[[x]])
  write.table(xx, paste0("data/fullsamples_test/Samples_iteration_gwr",sample(1:10000000, 1, replace=TRUE),".txt")) 
  
}

start <- Sys.time()
subsampling.plants(349113)
Sys.time() - start
# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns_raw <- fread("data/richness_patterns.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD, area) %>% 
  st_drop_geometry()

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

richness_patterns <- richness_patterns_raw %>% 
  left_join(continent_names, by = "LEVEL3_COD")


# Analysis ---------------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)



rich_overall_bru_bw <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by =c("LEVEL3_COD")) %>% 
  st_sf() %>% 
  as("Spatial") 

bw_global <- bw.gwr(sqrt(richness) ~ sqrt(richness2), data=rich_overall_bru_bw,
       approach = "AIC",
       adaptive = T)



rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL3_COD")) %>% 
  dplyr::select(-geometry)

setDT(plantlist_names)
setDT(dist_native)
setDT(rich_overall_bru)
setDT(midpoints_red)

output_model <- subsampling.plants(300000)
######


start <- Sys.time()
aa <- subsampling.plants(3000)
Sys.time() - start

plan(sequential)

Sys.time()
results <- future_lapply(seq(1,nrow(plantlist_names),1), subsampling.plants, future.seed = NULL, future.stdout = F)
Sys.time()
xx <- rbindlist(results)
write.table(xx, paste0("data/fullsamples_test/Samples_iteration_gwr_future_170K",sample(1:10000000, 1, replace=TRUE),".txt")) 


