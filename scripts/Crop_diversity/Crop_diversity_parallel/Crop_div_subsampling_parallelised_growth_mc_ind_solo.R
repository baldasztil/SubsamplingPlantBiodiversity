
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
options(dplyr.summarise.inform = FALSE)
# Defining functions -----------------------------------------------------------

# this is a function to subsample the WCVP
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
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
  
  dist <- dist_native_growth %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # richness patterns across tdwg3 areas
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3, growth_form) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = c("LEVEL_COD", "growth_form")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  # measuring correlation
  
  #correlation to overall richness 
  rich_rel_bru <- rich_rel %>% 
    ungroup() %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent correlation
  rich_rel_con <- rich_rel %>% 
    group_by(LEVEL1_COD) %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = as.character(LEVEL1_COD),
           sp = length(cumulative_namelist)) 
  
  if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
    print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  
  #extracting cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}

# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native_growth <-  fread("data/dist_native_growth.txt", fill = T) %>% 
  dplyr::select(-V1)
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 
growth <- fread("data/growth_rep_wcvp.txt", header = T, fill = T) %>% 
  dplyr::select(plant_name_id, growth_form)
# manipulate data --------------------------------------------------------------

# baseline richness ------------------------------------------------------------

# calculating overall richness patterns 

richness_patterns_con <- dist_native_growth %>% 
  group_by(continent_code_l1, growth_form) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- dist_native_growth %>% 
  group_by(region_code_l2, growth_form) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- dist_native_growth %>% 
  group_by(area_code_l3, growth_form) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns1$LEVEL_COD <- as.character(richness_patterns1$LEVEL_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns2$LEVEL_COD <- as.character(richness_patterns2$LEVEL_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)


# Analysis ---------------------------------------------------------------------


# creating objects for function
plantlist_dist <- dist_native_growth %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

######
for (i in 1:1)  {
  samples <- list()
  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(0,nrow(plantlist_names),1), mc.cores = 3, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/growth/Samples_iteration_",sample(1:10000000, 1, replace=TRUE),".txt")) # data/fullsamples_test/
}



