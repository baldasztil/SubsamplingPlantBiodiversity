library(tidyverse)
library(data.table)
library(sf)
library(parallel)

redlist_analysis <- read.csv("data/redlist/assessments.csv")
richness_patterns_allplants <- read.csv("output/Overall_richness.csv")

#Get the path of filenames
list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/redlist/SRLI/csv", full.names = TRUE)
#Read them in a list
samplelist <- lapply(list_data, fread)

redlist_index <- as.data.frame(do.call(rbind, samplelist))

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T)
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 

plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")


redlist_accepted <- plants_full_all %>% 
  filter(taxon_name %in% redlist_analysis$scientificName) %>% 
  left_join(redlist_analysis, by = c("taxon_name" = "scientificName")) 


redlist_plants_synonym <- wcvp_raw %>% 
  filter(taxon_name %in% redlist_analysis$scientificName) %>% 
  filter(!taxon_status =="Accepted") %>% 
  dplyr::select(plant_name_id, accepted_plant_name_id) %>% 
  mutate(accepted_plant_name_id = ifelse(accepted_plant_name_id =="",plant_name_id,accepted_plant_name_id)) %>% 
  filter(!accepted_plant_name_id %in% redlist_accepted$accepted_plant_name_id) 

redlist_plants_toadd <- plants_full_all %>% 
  filter(accepted_plant_name_id %in% redlist_plants_synonym$accepted_plant_name_id)

full_redlist <- bind_rows(redlist_accepted, redlist_plants_toadd)

write.csv(full_redlist, "redlist_full_syn.csv")



corrtax_srli_raw <- fread("output/redlist_srli/srli_noinfo_taxonomy_corr.csv",  header = T)
# cleaning 

corrtax_srli<- corrtax_srli_raw %>%
  mutate(taxon_name = str_trim(taxon_name))

corrtax <- wcvp_raw %>% 
  filter(taxon_name %in% corrtax_srli$taxon_name) 

hyb <- wcvp_raw %>% 
  filter(!species_hybrid == "") %>% 
  filter(taxon_name %in% corrtax_srli$taxon_name) 

corrtax <- wcvp_raw %>% 
  filter(taxon_name %in% corrtax_srli$taxon_name) 

corrtax_miss <- corrtax_srli %>% 
  filter(!taxon_name %in% wcvp_raw$taxon_name) 


srli_plants_accepted <- wcvp_raw %>% 
  filter(taxon_name %in% redlist_index$sp1) %>% 
  filter(taxon_status =="Accepted")

srli_plants_synonym_check <- wcvp_raw %>% 
  filter(taxon_name %in% redlist_index$sp1) %>% 
  filter(!taxon_status =="Accepted")

srli_plants_synonym <- wcvp_raw %>% 
  filter(taxon_name %in% redlist_index$sp1) %>% 
  filter(!taxon_status =="Accepted") %>% 
  dplyr::select(plant_name_id, accepted_plant_name_id) %>% 
  filter(!accepted_plant_name_id %in% srli_plants_accepted$accepted_plant_name_id) %>% 
  mutate(accepted_plant_name_id = ifelse(accepted_plant_name_id =="",plant_name_id,accepted_plant_name_id))
  

srli_plants_toadd <- plants_full_all %>% 
  filter(accepted_plant_name_id %in% srli_plants_synonym$accepted_plant_name_id)

srli_plants_corr_tax <- plants_full_all %>% 
  filter(taxon_name %in% corrtax_srli$taxon_name) %>% 
  filter(!taxon_name %in% srli_plants_toadd$taxon_name) %>% 
  filter(!taxon_name %in% srli_plants_accepted$taxon_name)

# all sp 
full_srli <- bind_rows(srli_plants_accepted, srli_plants_toadd, srli_plants_corr_tax)

write.csv(full_srli, "srli_full.csv")
dups <- full_srli %>% 
  group_by(plant_name_id) %>% 
  filter(n()>1)

# all accepted
srli_plants_accepted <- wcvp_raw %>% 
  filter(taxon_name %in% redlist_index$sp1) %>% 
  filter(taxon_status =="Accepted")

# all with no info on checklist 
srli_no_information_wcvp <- redlist_index %>% 
  filter(!sp1 %in% full_srli$taxon_name) %>%  
  filter(!sp1 %in% srli_plants_synonym_check$taxon_name) %>% 
  filter(!sp1 %in% corrtax_srli$sp1)

