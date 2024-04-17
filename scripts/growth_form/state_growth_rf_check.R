# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(ranger)
library(missRanger)
library(missForest)
library(mice)
library(miceRanger)
library(rcompanion)
library(vcd)
library(mltools)
library(regrrr)

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
gift_wood <- read.table("data/traits_all/woodiness_wcvp_data.txt", header = T, fill = T)  
gift_wood[gift_wood =="non-woody"] <- "herbaceous"


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
  filter(!location_doubtful == 1) %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)


dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            continent = paste(unique(continent), collapse = ','))

plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")



growths_rf <- read.csv("data/traits_all/growths_rf.csv") %>% 
  rename(certainty_growth = certainty)

state_rf <- read.csv("data/traits_all/woodi_herbaceous_rf.csv") %>% 
  rename(certainty_state = certainty)



full_data <- plants_full %>% 
  left_join(growths_rf, by = "plant_name_id") %>% 
  left_join(state_rf, by = "plant_name_id") 

tosave <- full_data %>% 
  dplyr::select(plant_name_id, family, genus, lifeform_description, state, certainty_state,  apd_plant_growth_form, certainty_growth)

write.csv(tosave, "data/traits_all/growth_state.csv")


checks <- full_data %>% 
  filter(genus == "Solanum") 

table(checks$state)
table(checks$apd_plant_growth_form)


checks2 <- full_data %>% 
  filter(family == "Araceae") %>% 
  filter(apd_plant_growth_form == "shrub")

state_check <- full_data %>% 
  left_join(gift_wood, by = "plant_name_id")

checks <- state_check %>% 
  filter(genus == "Viola") .

checks <- state_check %>% 
  filter(genus == "Solanum") %>% 
  filter(!lifeform_description == "")

table(checks$wood1)
table(checks$apd_plant_growth_form)
table(checks$lifeform_description)

state_check_check <- state_check %>% 
  filter(!lifeform_description == "")

xtabs(~ state + wood1, data = state_check_check)

b <-   as.data.frame(xtabs(~ state + wood1, data = state_check)) %>% 
  mutate(state = as.character(state), 
         wood1 = as.character(wood1)) %>% 
  mutate(true = state == wood1)  %>% 
  group_by(state) %>% 
  mutate(total = sum (Freq))
c <- b %>%
  filter(true == FALSE) %>%
  mutate(rel = round(Freq / total , 3))
#plot(b$rel ~ b$apd_plant_growth_form_origin)

d <-  c %>% 
  ungroup() %>% 
  group_by(state) %>% 
  reframe(rel_error = round(sum(Freq)/ total, 3)) %>% 
  filter(!duplicated(.))

e <- b %>% 
  ungroup() %>% 
  group_by(true) %>% 
  reframe(sum = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(total = sum(sum),
         rel_error = round(sum/total, 3))

