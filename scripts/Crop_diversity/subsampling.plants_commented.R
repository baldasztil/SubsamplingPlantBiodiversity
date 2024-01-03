
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

options(warn = -1) # turns of warnings if you want because at low species numbers you will get these
#options(warn = 0) # turns them back on
# function definition ----------------------------------------------------------
subsampling.plants <- function(spec_n) {
  
  #extract random sample of species from the overall species list
  species_sample <- sample_n(plantlist_names, spec_n)
  
  # this file contains the distribution of each species with one row corresponding
  # to one botanical country
  dist <- dist_native %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  #calculate richness patterns across botanical countries in sample
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3) %>% # group by botanical country (this can be ANY grouping)
    summarise(richness_sample = n_distinct(plant_name_id)) %>% # calculate species numbers
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0) %>% # replace NAs with 0 to include countries with no species
    mutate(sp = nrow(species_sample)) # indicator for species number
  
  # measuring correlation
  # calculate spearman and pearson coefficients 
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = as.numeric(cor(richness_sample, richness, 
                                      method="pearson"))) %>% 
    mutate(id = "overall",
           sp = nrow(species_sample)) 
  
  # does the same but for the patterns within each continent (again this 
  #can be ANY grouping)
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="pearson", data = ., exact =F)[[4]])  
  
  # merge the data sets 
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = nrow(species_sample)) %>% 
    left_join(spear, by = "id") 
  
  # you can use this to follow the progress when using lapply if you want
  #if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
  # print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  #}
  
  # cumulative pattern
  output_file <- bind_rows(rich_rel_bru, rich_rel_con) %>% 
    replace(is.na(.), 0) # replace NA correlation with 0 (no correlation)
  return(output_file)
}


# Import data ------------------------------------------------------------------

# this is the file with the species information
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
# this file contains the distribution of each species  
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 


# if you have species and distributions in one file it should also work! 
# Would just need to reformat some things / you could just generate a data frame
# that contains all the species names. 

# country codes 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

# shapefile with midpoints of all countries 
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)

# manipulate data --------------------------------------------------------------

# remove names without accepted name
wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

# remove all hybrids and genus names 
wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid =="")

# filter distrbutions 
dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% # not introduced 
  filter(!location_doubtful == 1) %>%  # not doubtful occurrence 
  filter(!area == "") %>% # not unknown
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)

# summarise information
dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

# create combined dataset
plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")

# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 

richness_patterns_con <- dist_native %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- dist_native %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)


# Analysis ---------------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native %>% 
  filter(plant_name_id %in% plantlist_names$plant_name_id) %>% 
  dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)


# global richness patterns as data frame 
rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>%  
  dplyr::select(-geometry)


# calculate correlations with a random sample at each step using lapply to do so sequentially
Sys.time()
samples <- lapply(seq(1,nrow(plantlist_names),1), subsampling.plants) 
Sys.time()

# export results as a .txt file 
xx <- rbindlist(samples)
write.table(xx, paste0("data/fullsamples_test/Samples_iteration_gwr_future_for",sample(1:10000000, 1, replace=TRUE),".txt")) 


