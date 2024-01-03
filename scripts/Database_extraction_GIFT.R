library(GIFT)
library(tidyverse)
library(data.table)
xx <- GIFT_traits_meta()

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 

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
  filter(!location_doubtful == 1) %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)


dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            continent = paste(unique(continent), collapse = ','))

plants_full_accepted <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")


woodiness <- GIFT_traits(trait_IDs = c("1.1.1"), agreement = 0.50000000000001,
                      bias_ref = FALSE, bias_deriv = FALSE)

#growth <- GIFT_traits(trait_IDs = c("1.2.1"), agreement = 0.50000001,
 #                    bias_ref = FALSE, bias_deriv = FALSE)

woodiness_dups <- woodiness %>% 
  filter(duplicated(work_species))

all <- plants_full_accepted %>%  
  filter(taxon_name %in% woodiness$work_species) 

notin <- woodiness %>%  
  filter(!work_species %in% plants_full_accepted$taxon_name) 


check <- wcvp_raw %>%  
  filter(taxon_name %in% notin$work_species) %>% 
  dplyr::select(taxon_name, accepted_plant_name_id) 

check_woodi <- notin %>% 
  left_join(check, by = c("work_species" = "taxon_name")) %>% 
  dplyr::select(trait_value_1.1.1, accepted_plant_name_id)

check_accepted <-plants_full_accepted %>% 
  filter(accepted_plant_name_id %in% check_woodi$accepted_plant_name_id) %>% 
  left_join(check_woodi, by = "accepted_plant_name_id") %>% 
  filter(!duplicated(.))


dupli <- check_accepted %>% 
  filter(duplicated(accepted_plant_name_id))


finetojoin <- check_accepted %>% 
  filter(!accepted_plant_name_id %in% dupli$accepted_plant_name_id) %>% 
  dplyr::select(accepted_plant_name_id, trait_value_1.1.1)


wood_fine <- woodiness %>% 
  filter(work_species %in% plants_full_accepted$taxon_name) %>% 
  dplyr::select(taxon_name = work_species, trait_value_1.1.1)


woodiness_full_raw <- plants_full_accepted %>% 
  left_join(finetojoin, by = "accepted_plant_name_id") %>% 
  left_join(wood_fine, by = c("taxon_name")) %>% 
  mutate(trait_value_1.1.1.x = ifelse(is.na(trait_value_1.1.1.x), trait_value_1.1.1.y,trait_value_1.1.1.x)) %>% 
  dplyr::select(wood1 = trait_value_1.1.1.x, plant_name_id)

write.csv(woodiness_full_raw, "data/traits_all/woodiness_gift_nofill.csv")

table(woodiness_full_raw$wood1)


growth_wcvp <- as.data.frame(unique(plants_full_accepted$lifeform_description))



plants_no_growth <- plants_full_accepted %>% 
  filter(plants_full_accepted$lifeform_description == "")

gen_list_wcvp <- plants_full_accepted %>% 
  group_by(family) %>% 
  reframe(genus = unique(genus))

woodiness_split <- woodiness %>% 
  mutate(taxon_name = work_species) %>% 
  separate(work_species, c("genus", "species"), sep = " ") 

woodiness_genus <- woodiness %>% 
  mutate(taxon_name = work_species) %>% 
  separate(work_species, c("genus", "species"), sep = " ") %>% 
  reframe(genus = unique(genus))


match <- woodiness_split %>% 
  filter(taxon_name %in% plants_no_growth$taxon_name)

nomatch <- plants_no_growth %>% 
  filter(!taxon_name %in% match$taxon_name)



fillable <- woodiness_genus %>% 
  filter(genus %in% plants_no_growth$genus)

genus_nogrowth <- gen_list_wcvp %>% 
  filter(genus %in% plants_no_growth$genus) %>% 
  filter(!genus %in% woodiness_genus$genus) 

non_fillable_sp <- genus_nogrowth %>% 
  filter(!genus %in% fillable$genus) %>% 
  inner_join(plants_no_growth, by ="genus",multiple = "all")


fillable_sp <- genus_nogrowth %>% 
  filter(genus %in% fillable$genus) %>% 
  inner_join(plants_no_growth, by ="genus",multiple = "all")


  

matching_wood <- wcvp_raw %>% 
  left_join(woodiness_split, by = "taxon_name") %>%
  group_by(accepted_plant_name_id) %>% 
  summarise(wood = paste(unique(trait_value_1.1.1), collapse = ','),
            plant_name_id = plant_name_id) %>% 
  separate(wood, c("wood1", "wood2", "wood3", "wood4"), sep = ",") 


wood_syns <- wcvp_raw %>% 
  left_join(woodiness_split, by = "taxon_name") %>%
  filter(taxon_status == "Synonym")

wood_syns_notpres <- plants_full_accepted %>% 
  filter(!plant_name_id %in% wood_syns$plant_name_id)


plant_wood <- plants_full_accepted %>% 
  left_join(matching_wood, by = "plant_name_id") %>% 
  mutate(wood1 = ifelse(wood1 == "NA", NA, wood1)) %>%
  mutate(wood2 = ifelse(wood2 == "NA", NA, wood2)) %>%
  mutate(wood1 = ifelse(is.na(wood1), wood2, wood1)) %>% 
  mutate(wood1 = ifelse(is.na(wood1), wood3, wood1)) %>% 
  mutate(wood1 = ifelse(is.na(wood1), wood4, wood1))

# THIS IS IMPORTANT ------------------------------------------------------------
# here you choose which species to fill information gaps for 

# all plants not in GIFT, but these might have information in the WCVP
plant_noinfo <- plant_wood %>% 
  filter(is.na(wood1))

# all plants were wood1 and 2 are not matching  
plant_notmaching <- plant_wood %>% 
  filter(!is.na(wood1) & !is.na(wood2)) %>% 
  filter(!wood1 == wood2)
           
plant_info <- as.data.frame(table(plant_wood$wood1))


#  all species no info on checklist or in GIFT
plant_checklistsub <- plants_full_accepted %>% 
  filter(plant_name_id %in% plant_noinfo$plant_name_id) %>% 
  filter(lifeform_description == "")

# ------------------------------------------------------------------------------
  

# freq of woody state per genus to fill species with no info
wood_freq <- plant_wood %>% 
  filter(genus %in% plant_noinfo$genus) %>%
  filter(!is.na(wood1)) %>% 
  group_by(genus) %>% 
  summarize(count_wood = length(which(wood1 == "woody")),
            count_nonwoody = length(which(wood1 == "non-woody")),
            count_variable = length(which(wood1 == "variable"))) 

# which one is the maximum of the columns 
gen_wood <- wood_freq %>% 
  mutate(max = case_when(max.col(wood_freq[2:4]) == 1 ~ 'woody', 
                         max.col(wood_freq[2:4]) == 2 ~ 'non-woody', 
                         max.col(wood_freq[2:4]) == 3 ~ 'variable'))



# genera with no info  
genus_nosource <- plant_checklistsub %>%
  filter(!genus %in% wood_freq$genus) %>% 
  reframe(genus = unique(genus))  
  

# species to fill in 
fill_sp <- plant_noinfo %>%  
  left_join(gen_wood, by = "genus") %>% 
  dplyr::select(plant_name_id, wood1 = max)

# ------------------------------------------------------------------------------

# filling in wood data for plants without info in checklist based on genus preference
wood_data <- plant_wood %>% 
  left_join(fill_sp, by = "plant_name_id") %>% 
  mutate(wood1.x = ifelse(is.na(wood1.x), wood1.y, wood1.x)) %>% 
  rename( wood1 = wood1.x) %>% 
  select(plant_name_id, wood1, wood2, wood3, wood4)


write.table(wood_data, "woodiness_wcvp_data.txt")
#table of scored spceis 
table(wood_data$wood1)

# lacking info overall 
plant_checklistsub <- plant_wood %>% 
  left_join(fill_sp, by = "plant_name_id") %>%  
  mutate(wood1.x = ifelse(is.na(wood1.x), wood1.y, wood1.x)) %>% 
  rename( wood1 = wood1.x) %>% 
  filter(is.na(wood1)) %>% 
  filter(lifeform_description == "")








wood_data_red <- wood_data %>% 
  dplyr::select(plant_name_id, wood1)
dist_native <- fread("data/dist_native.txt")

species_per_bru <- dist_native %>% 
  left_join(wood_data_red, by = "plant_name_id") %>% 
  group_by(area_code_l3) %>% 
  summarise(count_woodi = length(which(wood1 == "woody")),
            count_nonwoodi = length(which(wood1 == "non-woody")),
            sp = n_distinct(plant_name_id),
            diff = sp - (count_woodi + count_nonwoodi))



tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_3_names <- as.data.frame(tdwg_3) %>% 
  left_join(tdwg_1, by = "LEVEL1_COD") %>% 
  dplyr::select(LEVEL3_NAM, LEVEL3_COD, LEVEL1_COD, LEVEL1_NAM) %>% 
  summarise(continent_code_l1 = unique(LEVEL1_COD), 
            continent = unique(LEVEL1_NAM))



### woodi species 

  
species_per_cont <- dist_native %>% 
  group_by(continent) %>% 
  reframe(plant_name_id = unique(plant_name_id)) %>%
  left_join(wood_data_red, by = "plant_name_id") 
  
continent_summary <-   species_per_cont %>% 
  group_by(continent) %>% 
  summarise(count_woodi = length(which(wood1 == "woody")),
            count_nonwoodi = length(which(wood1 == "non-woody")),
            sp = n_distinct(plant_name_id),
            diff = sp - (count_woodi + count_nonwoodi),
            prop_nonwoodi = count_nonwoodi / (count_woodi + count_nonwoodi))
  
sum(species_per_cont$count_nonwoodi)  
  
length(which(is.na(species_per_cont$wood1) == T))
length(which(species_per_cont$wood1 == "variable"))
unique(species_per_cont$wood1)

species_per_cont 
