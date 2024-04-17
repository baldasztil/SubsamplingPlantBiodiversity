
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
library(ape)
library(phyloregion)
library(GWmodel)
library(feather)
library(vegan)
library(gmodels)
library(tmap)

dist_native <- dist_native <- fread("data/dist_native.txt") 
plantlist_names <- fread("data/wcvp_accepted_merged.txt")

gbif <- fread("data/dist_gbif.txt") 
srli <- read.csv("data/red/srli/srli_full.csv")
redlist <- read.csv("data/red/redlist_full_syn.csv")



dist_width <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(freq = n()) 

endemic <-  dist_width %>% 
  filter(freq == 1)

dist_width_gbif_wcvp <- gbif %>% 
  group_by(plant_name_id) %>% 
  summarise(freq_gbif = n()) %>% 
  right_join(dist_width, by = "plant_name_id") %>% 
  mutate(diff = freq - freq_gbif, 
         same = freq == freq_gbif) 

table(dist_width_gbif_wcvp$same)

dist_width_sum <- dist_width_gbif_wcvp %>%  
  group_by(diff) %>% 
  summarise(sp = n()) 

not_ingbif <- dist_native %>% 
  filter(!plant_name_id %in% gbif$plant_name_id)

aa <- gbif %>% 
  filter(!plant_name_id %in% dist_native$plant_name_id)


ranges_nogbif <- not_ingbif %>% 
  group_by(area_code_l3) %>% 
  summarise(n_miss = n_distinct(plant_name_id))

ranges_prop_missing <- dist_native %>% 
  group_by(area_code_l3) %>% 
  summarise(n = n_distinct(plant_name_id)) %>% 
  left_join(ranges_nogbif, by = "area_code_l3") %>% 
  mutate(diff = n_miss / n) %>% 
  left_join(tdwg_3, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  mutate_all(funs(replace_na(.,0))) %>% 
  st_as_sf()


map_missing_sp <- tm_shape(ranges_prop_missing) +
  tm_fill(col = "diff", palette = "viridis") +
  tm_borders(col = "black") +
  tm_layout(legend.outside = T, 
            title.snap.to.legend = F, 
            main.title = "Missing species in GBIF", 
            main.title.position = c("left", "top"),
            main.title.size = 1)  +
  tm_format("World")



dist_endemic <- dist_width_gbif_wcvp %>% 
  filter(freq == 1)

table(dist_endemic$diff)

dist_endemic_miss <- dist_width_gbif_wcvp %>% 
  filter(freq == 1 & is.na(freq_gbif))

endemics_nogbif <- not_ingbif %>% 
  filter(!plant_name_id %in% gbif$plant_name_id) %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  group_by(n) %>% 
  summarise(freq = n())

geog_endemics_miss <- dist_native %>% 
  filter(plant_name_id %in% dist_endemic_miss$plant_name_id) %>% 
  group_by(area_code_l3) %>% 
  summarise(number_missing = n_distinct(plant_name_id))


end_total <- dist_native %>% 
  filter(plant_name_id %in% dist_endemic$plant_name_id) %>% 
  group_by(area_code_l3) %>% 
  summarise(number = n_distinct(plant_name_id)) %>% 
  left_join(geog_endemics_miss, by = "area_code_l3") %>% 
  mutate(prop_miss = number_missing / number) %>% 
  right_join(tdwg_3, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  mutate_all(funs(replace_na(.,0))) %>% 
  st_as_sf()

map_missing_end_sp <- tm_shape(end_total) +
  tm_fill(col = "prop_miss", palette = "viridis") +
  tm_borders(col ="black") +
  tm_layout(legend.outside = T, 
            title.snap.to.legend = F, 
            main.title = "Missing endemics in GBIF", 
            main.title.position = c("left", "top"),
            main.title.size = 1)  +
  tm_format("World")
#map_missing_end_sp
missing_maps <-  tmap_arrange(map_missing_sp, map_missing_end_sp)

tmap_save(missing_maps, "wcvp_based_missing_sp.svg")

boxplot(end_total$prop_miss ~ end_total$HEMISPHERE)

  
geog_endemics_miss_cont <- dist_native %>% 
  filter(plant_name_id %in% dist_endemic_miss$plant_name_id) %>% 
  group_by(continent) %>% 
  summarise(number = n_distinct(plant_name_id))  %>% 
  mutate(number_tot = end_total$number, 
    diff = number / number_tot)

  
names_endemics_miss <- plants_full %>% 
  filter(plant_name_id %in% dist_endemic_miss$plant_name_id) %>% 
  group_by(genus) %>% 
  summarise(sp = n_distinct(plant_name_id))



ingbif <- dist_native %>% 
  filter(plant_name_id %in% gbif$plant_name_id)

ranges_ingbif <- ingbif %>% 
  group_by(plant_name_id) %>% 
  summarise(freq = n()) %>% 
  group_by(freq) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(sum = sum(n), 
         prop = n / sum)

ranges_total <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(freq = n()) %>% 
  group_by(freq) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  mutate(sum = sum(n), 
         prop = n / sum) %>% 
  left_join(ranges_ingbif, by = "freq") %>% 
  mutate(diff = round(prop.x - prop.y, 4))



head(ranges_total$freq, 10) - head(ranges_ingbif$freq, 10) 
