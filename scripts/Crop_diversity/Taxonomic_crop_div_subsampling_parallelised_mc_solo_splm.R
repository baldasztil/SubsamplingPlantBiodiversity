
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
  rich_rel <- dist %>% 
    group_by(LEVEL3_COD) %>% 
    summarise(richness_sample = n_distinct(family)) %>% 
    right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
    left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = nrow(species_sample),
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf() 
    
  

  global <- glance(spmodel::splm(richness ~ richness_sample, data=rich_rel,
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
    write.table(model_comp, paste0("data/taxonomic/checkpoints/Check@_sampling",
                                   nrow(species_sample),"sp.txt")) 
  }
  
  return(model_comp_out)
}
splm.continents <-  function(x, rich_rel)  {
  
  continent_data <- rich_rel %>% 
    filter(LEVEL1_NAM == x) 
  
  output_glance <- try(glance(splm(richness ~ richness_sample, data=continent_data,
                                   spcov_type =  c("exponential"))))
  
  if (class(output_glance)[1] =="try-error" & continent_data$sp[1] > 50000)  {
    output_glance_big <- try(glance(splm(richness ~ richness_sample, data=continent_data,
                                         spcov_type =  c("gravity"))))
    output_glance_big$continent <- x
    output_glance_big$sp <- unique(continent_data$sp)
    output_glance_big$model_id <- 2
    
    if (class(output_glance_big)[[1]] =="list")  {
      output_glance_fail <- global %>% 
        replace(!is.numeric(.) ==T, 1) %>% 
        mutate( 
          continent = unique(continent_data$LEVEL1_NAM),
          sp = unique(continent_data$sp),
          model_id = 0
        )
      output_glance_big <- output_glance_fail 
    }
      
    return(output_glance_big)
    
  }
  
  if (class(output_glance)[1] =="try-error" & continent_data$sp[1] < 50000)  {
    
    output_glance_fail <- global %>% 
      replace(!is.numeric(.) ==T, 0) %>% 
      mutate( 
        continent = unique(continent_data$LEVEL1_NAM),
        sp = unique(continent_data$sp),
        model_id = 0
      )
    return(output_glance_fail)
  }
  
  else {
    output_glance$continent <- x
    output_glance$sp <- unique(continent_data$sp)
    output_glance$model_id <- 1
    
  }
  
  return(output_glance)
}
# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full <- fread("data/wcvp_accepted_merged.txt")
richness_patterns_raw <- fread("data/richness_patterns_fam.txt")


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD) %>% 
  st_drop_geometry()

continent_names_vec <- unique(continent_names$LEVEL1_NAM)

richness_patterns <- richness_patterns_raw %>% 
  left_join(continent_names, by = "LEVEL3_COD")
  

# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)

plantlist_dist <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3, family)

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
subsampling.plants(349113)
Sys.time() - start

result_list <- lapply(seq(1,nrow(plantlist_names),1), 
                             subsampling.plants)
results <- rbindlist(result_list)
write.table(results, paste0("data/taxonomic/Samples_iteration_splm_",sample(1:10000000, 1, replace=TRUE),".txt")) 

