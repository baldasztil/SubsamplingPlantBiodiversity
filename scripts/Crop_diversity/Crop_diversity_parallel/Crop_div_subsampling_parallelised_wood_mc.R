
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
  
  dist <- dist_native %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # richness patterns across tdwg3 areas
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3, wood1) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = c("LEVEL_COD", "wood1")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  # measuring correlation
  
  #correlation to overall richness 
  rich_rel_bru <- rich_rel %>% 
    group_by(wood1) %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent correlation
  rich_rel_con <- rich_rel %>% 
    group_by(LEVEL1_COD, wood1) %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = as.character(LEVEL1_COD),
           sp = length(cumulative_namelist),
           wood1 = wood1) 
  
  if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
    print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  
  #extracting cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}

# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 
woodiness <- read.table("data/woodiness_wcvp_data.txt", header = T, fill = T)

# manipulate data --------------------------------------------------------------

wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

wcvp_sp <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(!taxon_rank %in% c("Genus")) 

wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid =="")

dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% 
  filter(!area == "") %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id) %>% 
  left_join(woodiness, by = "plant_name_id") %>% 
  filter(!is.na(wood1)) %>% 
  filter(!wood1 == "variable")

dist_nowood <- dist_native %>% 
  filter(is.na(wood1) | wood1 == "variable")

dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))




plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id") %>% 
  left_join(woodiness, by = "plant_name_id")  %>% 
  filter(!is.na(wood1)) %>% 
  filter(!wood1 == "variable")

plants_full_extinct_norange <- wcvp_accepted %>% 
  filter(!plant_name_id %in% plants_full$plant_name_id)

# baseline richness ------------------------------------------------------------

# calculating overall richness patterns 

richness_patterns_con <- dist_native %>% 
  group_by(continent_code_l1, wood1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2, wood1) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- dist_native %>% 
  group_by(area_code_l3, wood1) %>% 
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
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

###### running the function
for (i in 1:100)  {
  samples <- list()
  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(1,nrow(plantlist_names),1), mc.cores = 30, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/woodiness/Samples_iteration_wood_",i,".txt"))
}


