
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

options(dplyr.summarise.inform = FALSE)
# Defining functions -----------------------------------------------------------

# this is a function to subsample the WCVP

# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 
growth_form <- fread("data/traits_all/growth_forms.txt", header = T, fill = T)
wcvp_accepted <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")
# manipulate data --------------------------------------------------------------

name_match <- left_join(growth_form, wcvp_accepted, by = c("work_species" = "taxon_name"))



#na_groth <- is.na(growth_form)

growths <- growth_form %>% 
  filter(work_species %in% wcvp_accepted$taxon_name) 

inv_names <- growth_form %>% 
  filter(!work_species %in% wcvp_accepted$taxon_name)

sp_check <- wcvp_accepted %>% 
  filter(taxon_name %in% growth_form$work_species)

sp_nogrowth <- wcvp_accepted%>%  
  filter(!taxon_name %in% growth_form$work_species)

genus_nogrowth <- wcvp_accepted %>%  
  filter(genus %in% sp_nogrowth$genus) %>% 
  filter(!taxon_name %in% sp_nogrowth$genus)


wcvp_accepted_growth <- left_join(wcvp_accepted, growths, by = c("taxon_name" = "work_species"), multiple = "all")


# freq of woody state per genus to fill species with no info

wcvp_accepted_growth_rep <- data.frame(matrix(ncol = ncol(wcvp_accepted_growth), nrow = 0))
names(wcvp_accepted_growth_rep) <- names(wcvp_accepted_growth)

gen_names <- unique(wcvp_accepted_growth$genus)

# make progress bar for models 
n_iter <-  length(unique(wcvp_accepted_growth$genus)) # Number of iterations of the loop
# Initializes the progress bar
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = n_iter, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")  

for (i in 1 : length(unique(wcvp_accepted_growth$genus))) {
  setTxtProgressBar(pb,i)
  genus_accepted_growth <- wcvp_accepted_growth %>% 
    filter(genus == gen_names[i])
  
  gen_growth <- genus_accepted_growth %>% 
    filter(!is.na(growth_form))

  
  gen_nogrowth <- genus_accepted_growth %>% 
    filter(is.na(growth_form)) 
  
  if(nrow(gen_growth) == 0) {
    
    gen_nogrowth <- wcvp_accepted_growth %>% 
      filter(genus == gen_names[i])
    
    
    fam_growth <- wcvp_accepted_growth %>% 
      filter(family == gen_nogrowth$family[1]) %>% 
      filter(!is.na(growth_form))
    
    if (nrow(fam_growth) == 0) {
      rbind(wcvp_accepted_growth_rep, gen_nogrowth)
    }
    else {
      random_growth <- sample_n(fam_growth , nrow(gen_nogrowth), replace = T)
      
      gen_growth_rep <- gen_nogrowth %>% 
        mutate(growth_form = random_growth$growth_form)
      
      wcvp_accepted_growth_rep <- rbind(wcvp_accepted_growth_rep,gen_growth_rep) 
    }
  
  }
  else {
  random_growth <- sample_n(gen_growth , nrow(gen_nogrowth), replace = T)
  
  gen_growth_rep <- gen_nogrowth %>% 
    mutate(growth_form = random_growth$growth_form)
  
  full_gen <- rbind(gen_growth, gen_growth_rep)
  wcvp_accepted_growth_rep <- rbind(wcvp_accepted_growth_rep,full_gen)
  setTxtProgressBar(pb,i)
  }
}

check <- wcvp_accepted %>% 
  filter(plant_name_id %in% wcvp_accepted_growth_rep$plant_name_id)

miss <- wcvp_accepted_growth %>% 
  filter(!plant_name_id %in% check$plant_name_id)

miss$growth_form <- "subshrub"

full <- rbind(wcvp_accepted_growth_rep, miss)



growth_rest <- full %>% 
  distinct(plant_name_id, growth_form, .keep_all = TRUE) %>%  
  dplyr::select(plant_name_id, growth_form)

write.csv(growth_rest, "growth_rep_wcvp.txt") 

growth_rest <- read.csv("growth_rep_wcvp.txt")

dist_new <- dist_native %>% 
  left_join(growth_rest, by = "plant_name_id", multiple = "all")

write.csv(dist_new, "growth_distance.txt")


richness_patterns_con <- dist_new %>% 
  group_by(continent, growth_form) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent) %>% 
  mutate(ID = "con")

growth_rest <- read.csv("data/traits_all/growth_rep_wcvp.txt")
lifeform <- read.csv("data/traits_all/lifeform_mapping_original.csv")
wood <- fread("data/traits_all/woodiness_wcvp_data.txt")
growth_rest_comb <- growth_rest %>% 
  left_join(lifeform, by = "lifeform_description", multiple = "all") %>% 
  left_join(wood, by = "plant_name_id", multiple = "all")

flag <- growth_rest_comb %>% 
  filter(!is.na(wood2)) %>% 
  filter(!wood1 == wood2) %>% 
  mutate(flag = "flag") %>% 
  reframe(plant_name_id = plant_name_id, 
          flag = flag)

growth_cleaning <- growth_rest_comb %>% 
  left_join(flag, by = "plant_name_id", multiple = "all")

flag_cleaning <- growth_cleaning %>%
  filter(flag == "flag") %>% 
  group_by(lifeform_description, growth_form, humphreys_lifeform, wood1, wood2) %>% 
  reframe(n = n()) 

write.csv(flag_cleaning, "growth_form_flagcleaning.csv")

growth_inq <- growth_cleaning %>% 
  group_by(lifeform_description, growth_form, humphreys_lifeform, wood1, wood2) %>% 
  reframe(n = n()) 

growth_inq2 <- growth_rest_comb %>% 
  group_by(lifeform_description, growth_form, humphreys_lifeform, wood1, wood2) %>% 
  reframe(n = n()) 


write.csv(growth_inq2, "GIFT_WCVP_growth_form_unqiuecomb.csv")

growth_inq <- growth_rest %>% 
  group_by(lifeform_description, growth_form, life_form) %>% 
  reframe(n = n_distinct(.))


growth_inq <- growth_rest %>% 
  group_by(lifeform_description, growth_form, life_form) %>% 
  reframe(n = n_distinct(.))


