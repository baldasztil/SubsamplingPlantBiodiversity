
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
library(GWmodel)
library(feather)
library(ggsci)
library(plotrix)
library(rstatix)

# Defining functions -----------------------------------------------------------


subsampling.plants.mapping <- function(spec_n) {
  
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
           richness_sample = richness_sample,
           richness = richness) %>% 
    st_sf() 
  
  return(rich_rel)
}

# Import data ------------------------------------------------------------------
dist_native <- dist_native <- fread("data/dist_native.txt") 



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


# richness analysis
plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id, LEVEL3_COD = area_code_l3)
plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)

nrow(plantlist_names)

rich_overall_bru <- richness_patterns

rich_overall_bru_mid <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  mutate(richness2 = richness) %>% 
  left_join(midpoints_red, by = "LEVEL3_COD") %>% 
  st_as_sf()

rich_overall_bru_shp <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  dplyr::select(-LEVEL1_NAM) %>% 
  left_join(tdwg_3, by = "LEVEL3_COD") %>%
  mutate(richness2 = richness) %>% 
  st_as_sf() 

result_list_wcvp <- lapply(rep(2693,100), #nrow(plantlist_names)
                           subsampling.plants.mapping)

wcvp <- rbindlist(result_list_wcvp) %>% 
  mutate(dataset = "wcvp")


wcvp_stats <- wcvp %>%  
  group_by(LEVEL1_NAM) %>% 
  summarise(
    mean = mean(richness_sample),
    sd = sd(richness_sample),
    se = std.error(richness_sample), 
    
    median = median(richness_sample), 
    
    coef = sd / mean
    
  )






gbif_dist <- fread("data/dist_gbif.txt") 

plantlist_dist <- gbif_dist %>% 
  left_join(plants_full, by = c("plant_name_id")) %>% 
  dplyr::select(plant_name_id, LEVEL3_COD)



result_list_gbif <- lapply(rep(2693,100), #nrow(plantlist_names)
                           subsampling.plants.mapping)

gbif <- rbindlist(result_list_gbif) %>% 
  mutate(dataset = "gbif")

gbif_stats <- gbif %>%  
  group_by(LEVEL1_NAM) %>% 
  summarise(
    mean = mean(richness_sample),
    sd = sd(richness_sample),
    se = std.error(richness_sample), 
    
    median = median(richness_sample), 
    
    coef = sd / mean
    
  )
gbif_stats



data <- rbind(gbif, wcvp) %>% 
  filter(!LEVEL1_NAM == "ANTARCTICA") %>% 
  rename(continent = LEVEL1_NAM)

stats <- data %>%  
  group_by(LEVEL3_COD) %>% 
  wilcox_test(richness_sample ~ dataset, detailed = T)


ggplot(data, aes(x = reorder(continent,richness_sample), y =richness_sample, fill = dataset)) +
  #geom_jitter(color = "black", alpha = 0.05) +
  labs(x = "Continent", y = "richness") +
  scale_fill_viridis_d() +
  stat_boxplot(geom = "errorbar", width = 0.5, alpha = 0.95) +
  geom_boxplot(notch=F,alpha=0.95, width =0.5, color = "grey10", outlier.alpha = 0.05, outlier.size = 0.5) + 
  theme_bw()

stats_mapping <- rich_overall_bru_shp %>% 
  left_join(stats, by = "LEVEL3_COD") 
  
  
stats_mapping %>% 
  group_by(LEVEL1_NAM) %>% 
  summarise(median = mean(abs(estimate), na.rm = T))


ggplot(data = stats_mapping) +
  geom_sf(aes(fill = estimate)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick2", na.value = "white")  
  

plot(stats_mapping$estimate ~ stats_mapping$richness )
plot(stats_mapping$lat ~ stats_mapping$estimate )




rich_gbif <- plantlist_dist %>% 
  group_by(LEVEL3_COD) %>% 
  summarise(richness_sample = n_distinct(plant_name_id)) %>% 
  right_join(rich_overall_bru, by = "LEVEL3_COD") %>%
  left_join(midpoints_red, by = c("LEVEL3_COD")) %>% 
  replace(is.na(.), 0) %>% 
  mutate(sp = length(unique(plantlist_dist$plant_name_id)),
         richness_sample = richness_sample,
         richness = richness) %>% 
  st_sf() 

stats_gbif <- rich_gbif %>% 
  mutate(rich_diff = richness - richness_sample, 
         rich_diff_prop = 100- (richness_sample / richness *100) ) %>%
  mutate(rich_diff_prop = ifelse(is.infinite(rich_diff_prop), 100, rich_diff_prop)) %>% 
  st_drop_geometry()


stats_mapping_gbif <- rich_overall_bru_shp %>% 
  left_join(stats_gbif, by = "LEVEL3_COD") 



ggplot(data = stats_mapping_gbif) +
  geom_sf(aes(fill = rich_diff)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick2", na.value = "white")  


ggplot(data = stats_mapping_gbif) +
  geom_sf(aes(fill = rich_diff_prop)) +
  scale_fill_gradient(low = "lightgoldenrod", high = "firebrick2", na.value = "white")  

aa <- plantlist_dist %>% 
  filter(!paste(plant_name_id, LEVEL3_COD) %in% paste(dist_native$plant_name_id, dist_native$area_code_l3))
