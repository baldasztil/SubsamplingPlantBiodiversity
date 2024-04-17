
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(sf)
library(spmodel)



# Defining functions -----------------------------------------------------------
std.error <- function(x) sd(x)/sqrt(length(x))
splm.continents <-  function(x, rich_rel)  {
  
  continent_data <- rich_rel %>% 
    filter(LEVEL1_NAM == x) 
  
  output_glance <- try(glance(splm(richness ~ richness_sample, data=continent_data,
                               spcov_type =  c("exponential"))))
  
  if (class(output_glance)[1] =="try-error" & continent_data$sp[1] < 50000)  {
    output_glance <- global %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
             continent = unique(continent_data$LEVEL1_NAM),
             sp = unique(continent_data$sp),
             model_id = 0
      )
    return(output_glance)
  }
  
  if (class(output_glance)[1] =="try-error" & continent_data$sp[1] > 50000)  {
    output_glance <- try(glance(splm(richness ~ richness_sample, data=continent_data,
                                     spcov_type =  c("gravity"))))
    output_glance$continent <- x
    output_glance$sp <- unique(continent_data$sp)
    output_glance$model_id <- 2
    
    
    return(output_glance)
  }
  
  else {
    output_glance$continent <- x
    output_glance$sp <- unique(continent_data$sp)
    output_glance$model_id <- 1
    
  }
  return(output_glance)
}

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
  rich_rel <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample),
           richness_sample = richness_sample / area*1000000,
           richness = richness) %>% 
    st_sf() 
  hist(rich_rel$richness)
  global <- glance(splm(richness ~ richness_sample, data=rich_rel,
                      spcov_type =  c("exponential")))
  global$continent <- "GLOBAL"
  global$sp <- length(species_sample$plant_name_id)
  global$model_id <- 1
  
  continents <- rbindlist(lapply(continent_names_vec, splm.continents, 
                                 rich_rel = rich_rel ))
  
  model_comp <- rbind(global, continents)
  model_comp_out <- model_comp %>% 
    mutate_if(is.numeric,
                          round,
                          digits = 4)
  # cumulative pattern
  if (nrow(species_sample) %in% seq(1,nrow(plantlist_names),50000)){
    write.table(model_comp, paste0("data/richness/checkpoints/Check@_sampling",
                                   nrow(species_sample),"sp.txt")) 
  }
  
  return(model_comp_out)
}

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
  left_join(continent_names, by = "LEVEL3_COD") %>% 
  mutate(richness = richness / area*1000000)
  

# manipulate data --------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)

nrow(plantlist_names)

rich_overall_bru <- richness_patterns

rich_overall_bru_mid <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") 

rich_overall_bru_shp <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>% 
  st_as_sf()

setDT(plantlist_names)
setDT(dist_native)
setDT(rich_overall_bru)
setDT(midpoints_red)

global <- glance(splm(richness ~ 1, data = rich_overall_bru_shp, spcov_type = "exponential"))

# compare models ---------------------------------------------------------------

start <- Sys.time()
subsampling.plants(4000)
Sys.time() - start

write.table(aa, "test.txt")
result_list <- lapply(seq(1,nrow(plantlist_names),1), 
                             subsampling.plants)
results <- rbindlist(result_list)
write.table(results, paste0("data/richness/Samples_iteration_splm_",sample(1:10000000, 1, replace=TRUE),".txt")) 

