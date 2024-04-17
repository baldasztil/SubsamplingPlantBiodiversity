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
  left_join(wood_growth_gift_var, by = c("taxon_name" = "work_species")) %>% 
  mutate(wood1 = ifelse(is.na(wood1), state, wood1)) %>% 
  mutate(wood1 = ifelse(is.na(wood1), growth_gift, wood1), 
         lifeform_description = ifelse(lifeform_description == "", NA, lifeform_description))
  

xtabs(~state + wood1, wcvp_accepted_merger)

wcvp_accepted_merger_filt <- wcvp_accepted_merger %>% 
  filter(is.na(state))


checks <- wcvp_accepted_merger %>% 
  filter(state == "check")

checks <- wcvp_accepted_merger %>% 
  filter(is.na(lifeform_description))

checks <- wcvp_accepted_merger %>% 
  filter(is.na(wood1))


check_na_in(checks)

checks <- wcvp_accepted_merger %>% 
  filter(apd_plant_growth_form == "shrub" & state == "herbaceous") 


data <- wcvp_accepted_merger %>% 
  dplyr::select(lifeform_description, humphreys_lifeform, 
                wood1, state, apd_plant_growth_form, family, genus, plant_name_id, continent, growth_gift)

check_na_in(data)

data_ref <- data %>% 
  filter(state == "check" | is.na(state)) %>% 
  mutate(state = NA) 




data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) %>% 
  mutate(state = ifelse(apd_plant_growth_form == "shrub" & state == "herbaceous", "woody", state), 
         wood1 = ifelse(is.na(wood1), state, wood1)) %>% 
  mutate(wood1 = ifelse(is.na(wood1), "check", wood1))

data_merger_raw <- data_merger %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref)

check_na_in(data_merger_raw)

check_na_in(data_merger)

checks_lifeform_description <- data_ref %>% 
  filter(is.na(lifeform_description) & !is.na(wood1))

checks_wood1 <- data_ref %>% 
  filter(is.na(wood1) & !is.na(lifeform_description))
 
checks_both <- data_ref %>% 
  filter(is.na(wood1) & is.na(lifeform_description))

data_nona  <- data_merger %>% 
  filter(!is.na(state)) %>% 
  arrange(plant_name_id)

data_na_toinsert <- data_nona %>% 
  slice_sample(prop = nrow(data_ref) / nrow(data_merger)) %>% 
  arrange(plant_name_id)



na_lifeform_description <- round(nrow(data_na_toinsert) * (nrow(checks_lifeform_description) / nrow(data_ref)))
#na_wood1 <- round(nrow(data_na_toinsert) * (nrow(checks_wood1) / nrow(data_ref)))
na_both <- round(nrow(data_na_toinsert) * (nrow(checks_both) / nrow(data_ref)))



data_na <- data_na_toinsert %>%  
  mutate(state = NA)

data_na_growth <- data_na %>% 
  slice_sample(n = na_lifeform_description) %>% 
  mutate(lifeform_description = NA,
         humphreys_lifeform = NA)
  
data_na_both <- data_na %>% 
  filter(!plant_name_id %in% data_na_growth$plant_name_id) %>% 
  slice_sample(n = na_both) %>% 
  mutate(lifeform_description = NA,
         humphreys_lifeform = NA, 
    wood1 = NA)

data_na_full <- data_na %>% 
  filter(!plant_name_id %in% data_na_growth$plant_name_id & 
           !plant_name_id %in% data_na_both$plant_name_id) %>% 
  bind_rows(data_na_growth, data_na_both) %>% 
  arrange(plant_name_id)

check_na_in(data_ref)
check_na_in(data_na_full)

data_sim <- data_nona %>% 
  filter(!plant_name_id %in% data_na_full$plant_name_id) %>% 
  rbind(data_na_full) %>% 
  arrange(plant_name_id)

check_na_in(data_sim)

check_na_in(data_merger_raw)


data_origin <- data_nona


state <-  data_origin %>% 
  dplyr::select(state) 




data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, state, genus, continent, family) 

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)


mice_rf_imputed_v1 <- miceRanger(data_test_ref, 
                                 vars = c("state","wood1"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 1501
)



data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, state, genus,family, continent)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v2 <- miceRanger(data_test_ref, 
                                 vars = c("state"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 501
)






data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, state, genus, continent, family, wood1, lifeform_description)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v3 <- miceRanger(data_test_ref, 
                                 vars = c("state","wood1", "lifeform_description"),
                                 returnModels = T, 
                                 m = 1, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 500, 
                                 mtry = 2
)


nas <- mice_rf_imputed_v1$data %>% 
  filter(is.na(state))

data_test_na <- data_origin %>% 
  filter(plant_name_id %in% nas$plant_name_id) %>% 
  arrange(plant_name_id)



v1 <- as.data.frame(getVarImps(mice_rf_imputed_v1, 1:3, var = "state"))
v1$plant_name_id <- data_test_na$plant_name_id
v1$state_origin <- data_test_na$state
v1_imp <- v1 %>%  pivot_longer(cols = starts_with("Dataset"), values_to = "state_v1")

xtabs(~ Dataset_1 + state_origin, data = v1)
xtabs(~ Dataset_2 + state_origin, data = v1)
xtabs(~ Dataset_3 + state_origin, data = v1)


v2 <- as.data.frame(getVarImps(mice_rf_imputed_v2, 1:3, var = "state"))
v2_imp <- v2 %>%  pivot_longer(cols = starts_with("Dataset"), values_to = "state_v2", names_to = "name_v2")
v2$plant_name_id <- data_test_na$plant_name_id
v2$state_origin <- data_test_na$state

xtabs(~ Dataset_1 + state_origin, data = v2)
xtabs(~ Dataset_2 + state_origin, data = v2)
xtabs(~ Dataset_3 + state_origin, data = v2)


v3 <- as.data.frame(getVarImps(mice_rf_imputed_v3, 1:3, var = "state"))
v3_imp <- v3 %>%  pivot_longer(cols = starts_with("Dataset"), values_to = "state_v3", names_to = "name_v3")
v3$plant_name_id <- data_test_na$plant_name_id
v3$state_origin <- data_test_na$state

xtabs(~ Dataset_1 + state_origin, data = v3)
xtabs(~ Dataset_2 + state_origin, data = v3)
xtabs(~ Dataset_3 + state_origin, data = v3)



v4 <- as.data.frame(getVarImps(mice_rf_imputed_v4, 1:3, var = "state"))
v4_imp <- v4 %>%  pivot_longer(cols = starts_with("Dataset"), values_to = "state_v4", names_to = "name_v4")
v4$plant_name_id <- data_test_na$plant_name_id
v4$state_origin <- data_test_na$state

xtabs(~ Dataset_1 + state_origin, data = v4)
xtabs(~ Dataset_2 + state_origin, data = v4)
xtabs(~ Dataset_3 + state_origin, data = v4)


v5 <- as.data.frame(getVarImps(mice_rf_imputed_v5, 1:3, var = "state"))
v5_imp <- v5 %>%  pivot_longer(cols = starts_with("Dataset"), values_to = "state_v5", names_to = "name_v5")
v5$plant_name_id <- data_test_na$plant_name_id
v5$state_origin <- data_test_na$state

xtabs(~ Dataset_1 + state_origin, data = v5)
xtabs(~ Dataset_2 + state_origin, data = v5)
xtabs(~ Dataset_3 + state_origin, data = v5)





plotDistributions(mice_rf_imputed_v4)
plotCorrelations(mice_rf_imputed_v4)
plotVarConvergence(mice_rf_imputed_v4, vars = "allCategorical")
plotModelError(mice_rf_imputed_v4)
plotVarImportance(mice_rf_imputed_v4, display = "Relative")
plotImputationVariance(mice_rf_imputed_v4)

aa <- data_sim %>% 
  filter(genus == "Poa")



plotDistributions(mice_rf_imputed_v1)
plotCorrelations(mice_rf_imputed_v1)
plotVarConvergence(mice_rf_imputed_v1, vars = "allCategorical")
plotModelError(mice_rf_imputed_v1)
plotVarImportance(mice_rf_imputed_v1, display = "Relative")
plotImputationVariance(mice_rf_imputed_v1)






data_comp <- v4 %>% 
  left_join(data_test_na, by = "plant_name_id")
  

tab <- xtabs(~ state_origin + state_v1, data = v4)
tab <- xtabs(~ state_origin + state_v2, data = v4)
tab <- xtabs(~ state_origin + state_v3, data = v4)


summary(assocstats(tab))
mcc(v1_imp$state_v1, v2_imp$state_v2 ) # matthwees correlation coefficient

mice_rf_imputed_v1 <- mice_rf_imputed_v4
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


