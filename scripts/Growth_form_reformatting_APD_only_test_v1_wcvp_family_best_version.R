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



relative <- function(x) x/ rowSums(across(where(is.numeric)))

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
#woodiness <- read.table("data/traits_all/woodiness_wcvp_data.txt", header = T, fill = T)  

woodiness <- read.csv("data/traits_all/woodiness_gift_nofill.csv")
woodiness[woodiness =="non-woody"] <- "herbaceous"

lifeform <- read.csv("data/traits_all/lifeform_mapping_original.csv")

growthform <- fread("data/traits_all/traits_all.csv")



doubtful <- dist_raw %>% 
  filter(location_doubtful == 1)

tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
midpoints_red <- st_read(dsn ="data/wgsrpd-master/level3_midpoints") %>% 
  dplyr::select(LEVEL3_COD,geometry)

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

wcvp_accepted_merger <- plants_full %>%   
  left_join(lifeform, by = "lifeform_description") 


wcvp_accepted_merger_filt <- wcvp_accepted_merger %>% 
  filter(is.na(state))


check_na_in(wcvp_accepted_merger_filt)

checks <- wcvp_accepted_merger %>% 
  filter(apd_plant_growth_form == "shrub" & state == "herbaceous") 


data <- wcvp_accepted_merger %>% 
  dplyr::select(lifeform_description, 
                state, apd_plant_growth_form, family, genus, plant_name_id, continent)

check_na_in(data)

data_ref <- data %>% 
  filter(state == "check" | is.na(state)) %>% 
  mutate(state = NA, 
         apd_plant_growth_form = NA) 

data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) %>% 
  mutate(state = ifelse(apd_plant_growth_form == "shrub" & state == "herbaceous", "woody", state))




data_test <-  data_merger %>% 
  dplyr::select(plant_name_id, state, genus, continent, family) 

length(unique(data_test$plant_name_id))
table(data_test$state)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v1 <- miceRanger(data_test_ref, 
                                 vars = "state", 
                                 returnModels = T, 
                                 m = 99, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 501,
)



plotDistributions(mice_rf_imputed_v1)
plotCorrelations(mice_rf_imputed_v1)
plotVarConvergence(mice_rf_imputed_v1, vars = "allCategorical")
plotModelError(mice_rf_imputed_v1)
plotVarImportance(mice_rf_imputed_v1, display = "Relative")
plotImputationVariance(mice_rf_imputed_v1)


importance_frame <- data.frame()
for (i in 1: length(mice_rf_imputed_v1$allImps)) {
  data_list <- mice_rf_imputed_v1$allImport[[i]]
  data_frame <- do.call(rbind, data_list)
  data_frame$id <- paste0("dataset_",i)
  importance_frame <- rbind(importance_frame, data_frame)
  
}

importance_frame_rel <- importance_frame %>% 
  mutate(across(where(is.numeric), relative))


error_frame <- data.frame()
for (i in 1: length(mice_rf_imputed_v1$allError)) {
  data_frame <- mice_rf_imputed_v1$allError[[i]]
  data_frame$id <- paste0(i)
  error_frame <- rbind(error_frame, data_frame)
  
}


imps_frame <- data.frame()
for (i in 1: length(mice_rf_imputed_v1$finalImps)) {
  data_list <- mice_rf_imputed_v1$finalImps[i]
  names(data_list) <- "Dataset"
  data_frame <- as.data.frame(data_list$Dataset$state) %>% 
    rename(state = "data_list$Dataset$state")
  data_frame$id <- paste0(i)
  imps_frame <- rbind(imps_frame, data_frame)
}



best <- error_frame %>% 
  slice_max(state, n =1)


final_importance <- mice_rf_imputed_v1$finalImport
final_importance_rel <- do.call(rbind, final_importance) %>% 
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>%  
  mutate(across(where(is.numeric), relative))

#mice_rf_imputed_v1 <- paste0("mice_rf_imputed_v1$allImps$Dataset_",best$id,"$Iteration_",best$iteration)
#import <- rbindlist(mice_rf_imputed_v1$allError)
#imput <- as.data.frame(mice_rf_imputed_v1$finalImps$Dataset_1) # highest accuracy
#names(imput) <- "state_imp"


#table(imput$state_imp)

# imputed data
dataList <- completeData(mice_rf_imputed_v1)
data_imputed <- rbindlist(dataList)



nas <- data_test_ref %>% 
  filter(is.na(state))

imputed_statesv1 <- data_imputed %>% 
  filter(plant_name_id %in% nas$plant_name_id) %>% 
  dplyr::select(plant_name_id_v1 = plant_name_id, states_v1 = state)


length(unique(data_imputed$plant_name_id))
length(unique(data_test_ref$plant_name_id))

# checking for varying states
data_state <- data_imputed %>% 
  group_by(plant_name_id) %>% 
  mutate(count = length(unique(state)), 
         total = length(plant_name_id))

# filter multistate entries
data_multistate <- data_state %>% 
  filter(count == 2)

length(unique(data_multistate$plant_name_id))

# filter unistate entries and remove geographic duplicates
data_unistate_small <- data_state %>% 
  filter(count == 1) 

length(unique(data_unistate_small$plant_name_id))

# look into multistate data to determine which ones to keep
data_multistate_exam <- data_state %>%
  ungroup() %>% 
  group_by(plant_name_id, state) %>% 
  mutate(freq = length(plant_name_id), 
         rel = freq / total) %>% 
  ungroup() %>% 
  filter(!duplicated(.))

# look at the once with more than 50% in one category
data_multistate_rel_adj <- data_multistate_exam %>% 
  group_by(plant_name_id) %>% 
  slice_max(rel, n = 1) %>% 
  filter(rel > 0.5) %>% 
  dplyr::select(plant_name_id, state, certainty = rel)

save.image("99_datasets.RData")
write.csv(data_multistate_rel_adj, "data/traits_all/woodi_herbaceous_rf.csv")

