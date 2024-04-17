
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(vegan)
library(picante)


# Import data ------------------------------------------------------------------
# this is the file with the species information
wcvp_raw <- fread("data/wcvp/wcvp_names_032023.csv", 
                  sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8") 
# this file contains the distribution of each species  
dist_raw <- fread("data/wcvp/wcvp_distribution_032023.csv", 
                  sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")


midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)


# geographic data --------------------------------------------------------------
tdwg_3 <- rWCVP::wgsrpd3
tdwg_3$area <- as.numeric(st_area(tdwg_3_full) / 1000000)


tdwg_3_mapping <- rWCVP::wgsrpd_mapping %>% 
  filter(!duplicated(LEVEL3_COD)) %>% 
  dplyr::select(HEMISPHERE, LEVEL1_NAM, LEVEL2_NAM, LEVEL3_COD, COUNTRY)

tdwg_3_full <- tdwg_3 %>% 
  left_join(tdwg_3_mapping, by = "LEVEL3_COD")


midpoints_raw <- st_centroid(tdwg_3)

midpoints_red <- midpoints_raw %>% 
  dplyr::select(LEVEL3_COD,geometry)

tdwg_3_midpoints <- midpoints_raw %>% 
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  dplyr::select(LEVEL3_COD, lon,lat) %>% 
  left_join(tdwg_3_full, by = "LEVEL3_COD" )

write_sf(tdwg_3_midpoints, "data/wgsrpd-master/level3/tdwg_3.shp")
write_sf(midpoints_red, "data/wgsrpd-master/level3_midpoints/midpoints.shp")

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3")
midpoints <- st_read(dsn ="data/wgsrpd-master/level3_midpoints")


# baseline patterns ------------------------------------------------------------

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


dist_1 <- dist_raw %>% 
  filter(introduced == 0) %>% 
  filter(!location_doubtful == 1) %>% 
  filter(!extinct == 1) %>% 
  filter(!area == "") %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)
  

dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))


plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")


dist_native <- dist_1 %>% filter(plant_name_id %in% plants_full$plant_name_id)
length(unique(dist_native$plant_name_id))
write.table(dist_native, "data/dist_native.txt")
dist_native <- fread("data/dist_native.txt")

length(unique(dist_1$plant_name_id))

wcvp_accepted_merged <- plants_full %>% filter(plant_name_id %in% dist_native$plant_name_id)
length(unique(wcvp_accepted_merged$plant_name_id))
write.table(wcvp_accepted_merged, "data/wcvp_accepted_merged.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")

#plants_full_extinct_norange <- wcvp_accepted %>% 
 # filter(!plant_name_id %in% plants_full$plant_name_id)

# baseline richness ------------------------------------------------------------

dist_native_calc <- dist_native %>% 
  left_join(plants_full, by = "plant_name_id")

richness_patterns_bru <- dist_native_calc %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns_bru_fam_rep <- dist_native_calc %>% 
  group_by(area_code_l3, family) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")


write.table(richness_patterns_bru_fam_rep, "richness_patterns_fam.txt")

write.table(richness_patterns_bru, "richness_patterns.txt")


rich_overall_bru <- richness_patterns_bru %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL3_COD")) %>% 
  dplyr::select(-geometry)

# baseline patterns ------------------------------------------------------------

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(midpoints_red, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  dplyr::select(-geometry)