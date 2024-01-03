
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------

# This is a script to create random subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species richness in the subsample to global patterns.

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(parallel)
library(doParallel)
# Defining functions -----------------------------------------------------------

# this is a function to subsample the WCVP
extracting.plants <- function(spec_n) {
  list_rich_cumulative <- list()
  list_rel_cumulative <- list()
  
  for (i in 1:100) { 
  print(paste0("This is iteration ", i))
  cumulative_namelist <- c()
  # accumulative subsampling 
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
  
  # richness patterns across tdwg3 areas
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  # measuring correlation
  
  #correlation to overall richness 
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent correlation
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="pearson", data = ., exact =F)[[4]])  
  
  # combining the results
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  
  #extracting cumulative pattern
  list_rel_cumulative[[i]] <- bind_rows(rich_rel_bru, rich_rel_con)
  list_rich_cumulative[[i]] <- rich_rel
  }
  assign("list_rel_cumulative", list_rel_cumulative, envir = .GlobalEnv)
  return(list_rich_cumulative)
}

# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------
richness_patterns <- read.csv("output/Overall_richness.csv")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

# Analysis ---------------------------------------------------------------------


# creating objects for function
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))


xx <- extracting.plants(1390)

###### running the function
for (i in 1:100)  {
  samples <- list()
  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(0,nrow(plantlist_names),1), mc.cores = 30, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/fullsamples/Samples_iteration_end_",i,".txt"))
}


