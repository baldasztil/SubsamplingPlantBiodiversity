
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
library(data.table)



# Defining functions -----------------------------------------------------------
std.error <- function(x) sd(x)/sqrt(length(x))
splm.continents <-  function(x, rich_rel)  {
  
  continent_data <- rich_rel %>% 
    filter(LEVEL1_NAM == x) 
  
  output_glance <- try(glance(splm(richness ~ richness_sample, data=continent_data,
                               spcov_type =  c("gravity"))))
  
  if (class(output_glance)[1] =="try-error")  {
    output_glance <- global %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate(sp = unique(continent_data$sp), 
             continent = unique(continent_data$LEVEL1_NAM))
    return(output_glance)
    }
  else {
    output_glance$continent <- x
    output_glance$sp <- unique(continent_data$sp)
  }
  return(output_glance)
}

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
  
  dist <- plantlist_dist %>% 
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
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf() 
    
    rich_rel_shp <- rich_rel %>% 
    as("Spatial")
  

  global <- glance(spmodel::splm(richness ~ richness_sample, data=rich_rel,
                      spcov_type =  c("gravity")))
  global$continent <- "GLOBAL"
  global$sp <- length(species$plant_name_id)
  
  
  
  continents <- rbindlist(lapply(continent_names_vec, splm.continents, 
                                 rich_rel = rich_rel ))
  
  model_comp <- rbind(global, continents)
  model_comp$sp <- length(species$plant_name_id)
  
  
  return(model_comp)
}









# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") %>% 
  rename(area_)
plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns_raw <- fread("data/richness_patterns.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

richness_patterns <- richness_patterns_raw %>% 
  left_join(continent_names, by = "LEVEL3_COD")
  

# manipulate data --------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)

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

start <- Sys.time()
result_list <- future_lapply(seq(1,10000,10), 
                             subsampling.plants, 
                             future.seed = NULL, future.stdout = F)
Sys.time() - start

aa <- rbindlist(result_list)
bb <- rbindlist(result_list)
aa - bb



setDT(plantlist_names)
setDT(dist_native)
setDT(rich_overall_bru)
setDT(midpoints_red)
subsampling.plants.dt <- function(spec_n) {
  cumulative_namelist <- c()
  
  # Convert data frames to data.tables
  
  while (TRUE) {
    plantlist_names_left <- plantlist_names[!plant_name_id %in% cumulative_namelist]
    
    if (nrow(plantlist_names_left) > spec_n) {
      species_sample <- plantlist_names_left[sample(.N, spec_n)]
    } else {
      species_sample <- plantlist_names_left
    }
    rich_rel <- dist[, .(richness_sample = uniqueN(plant_name_id)), by = .(LEVEL3_COD)]
    
    cumulative_namelist <- c(cumulative_namelist, species_sample$plant_name_id)
    
    species <- plantlist_names[plant_name_id %in% cumulative_namelist]
    
    dist <- plantlist_dist[plant_name_id %in% species$plant_name_id]
    
    rich_rel <- dist[, .(richness_sample = uniqueN(plant_name_id)), by = LEVEL3_COD][
      rich_overall_bru, on = "LEVEL3_COD"][midpoints_red, on = "LEVEL3_COD"][, .(LEVEL3_COD, LEVEL1_NAM,  richness_sample, richness, geometry)][is.na(richness_sample), 
       richness_sample := 0][, sp := length(cumulative_namelist)]
    
    rich_rel <- rich_rel %>% 
      st_as_sf() 
    

    global <- glance(spmodel::splm(richness ~ richness_sample, data = rich_rel, 
                                   spcov_type = c("gravity")))
    global$continent <- "GLOBAL"
    global$sp <- nrow(species)
    
    continents <- rbindlist(lapply(continent_names_vec, splm.continents, rich_rel = rich_rel))
    
    model_comp <- rbindlist(list(global, continents))
    model_comp$sp <- nrow(species)
    
    return(model_comp)
  }
}

library(rbenchmark)

benchmark(
  "dt" = {
 x <- lapply(sample(1:length(plantlist_names$plant_name_id), 5), subsampling.plants.dt)
},
 "nodt" = {
  x <- lapply(sample(1:length(plantlist_names$plant_name_id), 5), subsampling.plants)
},
replications = 100,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))

length(plantlist_names$plant_name_id)
start <- Sys.time()
result_list <- lapply(sample(1:length(plantlist_names$plant_name_id), 5), subsampling.plants)
Sys.time() - start
