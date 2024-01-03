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




relative <- function(x) x/ rowSums(across(where(is.numeric)))

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
woodiness <- read.table("data/traits_all/woodiness_wcvp_data.txt", header = T, fill = T)  
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

wood_growth_gift <- growthform %>%  
  mutate(growth_gift = case_when(growth_form %in% c("subshrub", "shrub","tree") ~ "woody", 
                               growth_form %in% c("herb") ~ "herbaceous", 
                               TRUE ~ NA
         )) %>% 
  dplyr::select(growth_gift, work_species) %>% 
  filter(!duplicated(.))


names_dup <- wood_growth_gift %>% 
  filter(!is.na(growth_gift)) %>% 
  group_by(work_species) %>% 
  summarise(unique = length(unique(growth_gift))) %>% 
  filter(unique > 1)


wood_growth_gift_var <- wood_growth_gift %>%  
  mutate(growth_gift = ifelse(work_species %in% names_dup$work_species, "variable", growth_gift)) %>% 
  dplyr::select(growth_gift, work_species) %>% 
  filter(!duplicated(.)) %>% 
  filter(!is.na(growth_gift)) %>% 
  group_by(work_species) %>% 
  mutate(count = length(work_species))

length(unique(wood_growth_gift_var$work_species))

unique(wood_growth_gift_var$growth_gift)

wcvp_accepted_merger <- plants_full %>%   
  left_join(lifeform, by = "lifeform_description") %>% 
  left_join(woodiness, by = "plant_name_id") %>% 
  left_join(wood_growth_gift_var, by = c("taxon_name" = "work_species"))

wcvp_accepted_merger_filt <- wcvp_accepted_merger %>% 
  filter(is.na(state))


checks <- wcvp_accepted_merger %>% 
  filter(apd_plant_growth_form == "check")

checks <- wcvp_accepted_merger %>% 
  filter(lifeform_description == "")

table(checks$growth_gift)
table(checks$wood1)

xtabs(~ wood1 + growth_gift, data = checks)


checks <- wcvp_accepted_merger %>% 
  filter(apd_plant_growth_form == "shrub" & state == "herbaceous") 


data <- wcvp_accepted_merger %>% 
  dplyr::select(lifeform_description, humphreys_lifeform, 
                wood1, state, apd_plant_growth_form, family, genus, plant_name_id, continent, growth_gift)


data_ref <- data %>% 
  filter(state == "check" | is.na(state)) %>% 
  mutate(state = NA) 


data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) %>% 
  mutate(humphreys_lifeform = ifelse(is.na(humphreys_lifeform), "", humphreys_lifeform), 
         state = ifelse(apd_plant_growth_form == "shrub" & state == "herbaceous", "woody", state), 
         wood1 = ifelse(is.na(wood1), state, wood1)) %>% 
  mutate(wood1 = ifelse(is.na(wood1), "variable", wood1))


data_nona <- data_merger %>% 
  mutate(apd_plant_growth_form = ifelse(is.na(state),NA, apd_plant_growth_form), 
         lifeform_description = ifelse(is.na(state),NA, lifeform_description))



data_test <-  data_nona %>% 
  dplyr::select(plant_name_id, state, genus, continent, family, growth_gift, wood1) %>% 
  mutate(growth_gift = ifelse(is.na(growth_gift),"", growth_gift), 
         wood1 = ifelse(is.na(wood1),"", wood1))

length(unique(data_test$plant_name_id))
xtabs(~ growth_gift + wood1, data_test)


data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v1 <- miceRanger(data_test_ref, 
                                 vars = "state", 
                                 returnModels = T, 
                                 m = 10, 
                                 maxiter = 10, 
                                 valueSelector = "value", 
                                 num.trees = 1000,
)




names <- mice_rf_imputed_v1$data$plant_name_id
nas <-as.data.frame(mice_rf_imputed_v1$naWhere)$state


data_test_na <- data_origin %>% 
  cbind(nas) %>% 
  filter(nas == TRUE) 

v1 <- as.data.frame(getVarImps(mice_rf_imputed_v1, 1:3, var = "state"))
v1$plant_name_id <- data_test_na$plant_name_id
v1$state_origin <- data_test_na$state
v1_imp <- v1 %>%  pivot_longer(cols = starts_with("Dataset"), values_to = "state_v1")

xtabs(~ Dataset_1 + state_origin, data = v1)
xtabs(~ Dataset_2 + state_origin, data = v1)
xtabs(~ Dataset_3 + state_origin, data = v1)



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
dataListv2 <- completeData(mice_rf_imputed_v2)
data_imputedv2 <- rbindlist(dataListv2)



nas <- data_test_ref %>% 
  filter(is.na(state))

imputed_statesv1 <- data_imputed %>% 
  filter(plant_name_id %in% nas$plant_name_id) %>% 
  dplyr::select(plant_name_id_v1 = plant_name_id, states_v1 = state)

imputed_statesv2 <- data_imputedv2 %>% 
  filter(plant_name_id %in% nas$plant_name_id) %>% 
  dplyr::select(plant_name_id_v2 = plant_name_id, states_v2 = state)

imputed_states_comp <- cbind(imputed_statesv1, imputed_statesv2)
imputed_states_comp$name_id <- imputed_states_comp$plant_name_id_v1 ==imputed_states_comp$plant_name_id_v2

tab <- xtabs(~states_v1 + states_v2 , data = imputed_states_comp)


length(unique(data_imputed$plant_name_id))
length(unique(data_test_ref$plant_name_id))

# checking for varying states
data_state <- data_imputed %>% 
  group_by(plant_name_id) %>% 
  mutate(count = length(unique(state))) %>%  
  mutate(total = length(plant_name_id)) 

# remove geographic duplicates 
data_state_small <- data_state %>% 
  ungroup() %>%  
  dplyr::select(plant_name_id, state, count) %>% 
  filter(!duplicated(.))

length(unique(data_state_small$plant_name_id))

# filter multistate entries
data_multistate <- data_state_small %>% 
  filter(count == 2)

length(unique(data_multistate$plant_name_id))

# filter unistate entries and remove geographic duplicates
data_unistate_small <- data_state_small %>% 
  filter(count == 1) 

length(unique(data_unistate_small$plant_name_id))

# look into multistate data to determine which ones to keep
data_multistate_exam <- data_state %>%
  ungroup() %>% 
  group_by(plant_name_id) %>% 
  mutate(total = length(plant_name_id)) %>% 
  group_by(plant_name_id, state) %>% 
  mutate(freq = length(plant_name_id)) %>% 
  group_by(plant_name_id) %>% 
  mutate(rel = freq/total)


data_tables <- data_multistate_exam %>%
  ungroup() %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, state) %>% 
  filter(!duplicated(.))

tab <- xtabs(~ state + apd_plant_growth_form , data = data_tables)

data_impt_exam <- data_multistate_exam %>% 
  ungroup() %>% 
  dplyr::select(-continent) %>% 
  dplyr::select(-apd_plant_growth_form) %>% 
  filter(!duplicated(.))


# look at the once with more than 50% in one category
data_multistate_rel <- data_impt_exam %>% 
  group_by(plant_name_id) %>% 
  slice_max(rel, n = 1) %>% 
  filter(rel > 0.5)



# ectract all 
data_multistate_exam_split <- data_impt_exam %>% 
  filter(rel== 0.5)

huh <- data_multistate_exam_split %>% 
  filter(duplicated(plant_name_id))

huh2 <- data_multistate_exam_split %>% 
  filter(!duplicated(plant_name_id))

huh3 <- huh %>% 
  filter(duplicated(plant_name_id))

length(unique(data_multistate_exam_split$plant_name_id))

# association analysis of the two variables
data_na_test <- data_test_ref %>% 
  filter(is.na(state))

data_origin_stat <- data_test %>% 
  filter(plant_name_id %in% data_na_test$plant_name_id) 

data_state_small_real <- data_origin_stat %>% 
  dplyr::select(plant_name_id, state) %>% 
  filter(!duplicated(plant_name_id))

cramerxx <- data_imputed %>% 
  ungroup() %>%  
  dplyr::select(plant_name_id, continent, imp_state = state) %>% 
  filter(!duplicated(.)) %>% 
  right_join(data_state_small_real, by = "plant_name_id")


tab <- xtabs(~state + imp_state , data = cramerxx)

summary(assocstats(tab))
mosaicplot(tab)

library(mltools)

actuals <- as.factor(cramerxx$state)
preds <- as.factor(cramerxx$imp_state)
mcc(preds, actuals) # matthwees correlation coefficient


