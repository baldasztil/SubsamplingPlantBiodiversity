library(tidyverse)
library(data.table)

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
plants_full <- fread("data/wcvp_accepted_merged.txt")
plants_full_accepted <- fread("data/wcvp_accepted_merged.txt")


#biome <- GIFT_traits(trait_IDs = c("1.2.1"), agreement = 0.51,
                   #   bias_ref = FALSE, bias_deriv = FALSE)


biome_wcvp <- as.data.frame(unique(plants_full_accepted$climate_description))

plants_no_biome <- plants_full_accepted %>% 
  filter(plants_full_accepted$climate_description == "")

plants_biome<- plants_full_accepted %>% 
  filter(!plants_full_accepted$climate_description == "") %>%  
  mutate(climate_description = ifelse(climate_description == "Desert or Dry Shrubland", "Desert and Dry Shrubland", climate_description)) 



gen_list_wcvp <- plants_full_accepted %>% 
  group_by(family) %>% 
  reframe(genus = unique(genus))


gen_list_nobiome <- plants_no_biome %>% 
  group_by(family) %>% 
  reframe(genus = unique(genus))



plants_gen_nobiome <- plants_full_accepted %>% 
  filter(genus %in% gen_list_nobiome$genus)

# THIS IS IMPORTANT ------------------------------------------------------------
# here you choose which species to fill information gaps for 

# ------------------------------------------------------------------------------
plants_replace <- plants_gen_nobiome %>% 
mutate(climate_description = ifelse(climate_description == "Desert or Dry Shrubland", "Desert and Dry Shrubland", climate_description)) 

unique(plants_replace$climate_description)


plants_replace_filt <- plants_full_accepted %>% 
  filter(genus == "Acritopappus")

# freq of woody state per genus to fill species with no info
biome_freq_climate <- plants_replace %>% 
  group_by(genus) %>% 
  summarize(Seasonally_Dry_Tropical = length(which(climate_description == "Seasonally Dry Tropical")),
            Wet_Tropical = length(which(climate_description == "Wet Tropical")),
            Desert_and_Dry_Shrubland = length(which(climate_description == "Desert and Dry Shrubland")),
            Subtropical = length(which(climate_description == "Subtropical")),
            Temperate = length(which(climate_description == "Temperate")),
            Montane_Tropical = length(which(climate_description == "Montane Tropical")),
            Subalpine_or_Subarctic = length(which(climate_description == "Subalpine or Subarctic")), 
            Subtropical_and_Tropical = length(which(climate_description == "Subtropical and Tropical"))) 


#sample(c("A", "B", "C"), nrow(df), prob = c(0.2, 0.3, 0.5), replace = TRUE)


biome_freq <- biome_freq_climate %>% 
  mutate(sum = rowSums(biome_freq_climate[2:9]))



# which one is the maximum of the columns 

gen_biome <- biome_freq %>% 
  mutate(max = colnames(biome_freq[2:9])[max.col(biome_freq[2:9], ties.method = "random")]) 




# species to fill in 
fill_sp <- plants_no_biome %>%  
  left_join(gen_biome, by = "genus") 


fill_sp_red <- fill_sp %>% 
  filter(sum > 0) %>% 
  dplyr::select(plant_name_id, climate_description_fill = max) %>% 
  mutate(climate_description_fill = gsub("_", " ", climate_description_fill))


fill_sp_na <- fill_sp %>% 
  filter(sum == 0)

table <- plants_biome$climate_description

for (i in 1:nrow(fill_sp_na)) {
  fill_sp_na$climate_description_fill[i] <- sample(table, 1, replace = F)
}



fill_sp_full <- bind_rows(fill_sp_red,fill_sp_na) %>% 
  dplyr::select(plant_name_id, climate_description_fill)

# ------------------------------------------------------------------------------

# filling in wood data for plants without info in checklist based on genus preference
biome_data <- plants_full_accepted %>% 
  left_join(fill_sp_full, by = "plant_name_id") %>% 
  mutate(climate_description = ifelse(climate_description == "", climate_description_fill, climate_description)) 

biome_data_red <- biome_data %>% 
  dplyr::select(plant_name_id, biome = climate_description)

write.table(biome_data_red, "biome_wcvp_data.txt")


#table of scored species 
table(biome_data$climate_description)

length(which(is.na(biome_data$climate_description) == FALSE))

length(unique(biome_data$plant_name_id))



dist_native <- fread("data/dist_native.txt")

species_per_bru <- dist_native %>% 
  left_join(biome_data_red, by = "plant_name_id") %>% 
  reframe(sp_per_biome = table(climate_description), .by = c("area_code_l3", "climate_description"))


test <- species_per_bru %>% 
  filter(area_code_l3 == "AUT")

library(sf)
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_3_names <- as.data.frame(tdwg_3) %>% 
  left_join(tdwg_1, by = "LEVEL1_COD") %>% 
  dplyr::select(LEVEL3_NAM, LEVEL3_COD, LEVEL1_COD, LEVEL1_NAM) %>% 
  summarise(continent_code_l1 = unique(LEVEL1_COD), 
            continent = unique(LEVEL1_NAM))


tdwg_3_names_j <- tdwg_3_names %>% 
  mutate(continent_code_l1 = as.character(continent_code_l1))

richness <- richness_patterns %>% 
  right_join(tdwg_3_names_j, by = c("LEVEL_COD" = "continent_code_l1"))

### woodi species 

  
species_per_cont <- dist_native %>% 
  group_by(continent) %>% 
  reframe(plant_name_id = unique(plant_name_id)) %>%
  left_join(biome_data_red, by = "plant_name_id") 
  
continent_summary <-   species_per_cont %>% 
  reframe(sp_per_biome = table(climate_description), .by = c("continent", "climate_description"))

