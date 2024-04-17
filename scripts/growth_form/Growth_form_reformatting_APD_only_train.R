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



relative <- function(x) x/ rowSums(across(where(is.numeric)))

# Import data ------------------------------------------------------------------
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
woodiness <- read.table("data/traits_all/woodiness_wcvp_data.txt", header = T, fill = T)  
woodiness[woodiness =="non-woody"] <- "herbaceous"

lifeform <- read.csv("data/traits_all/lifeform_mapping_original.csv")

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
 

wcvp_accepted_merger <- wcvp_accepted %>%   
  left_join(lifeform, by = "lifeform_description") %>% 
  left_join(woodiness, by = "plant_name_id")

wcvp_accepted_merger_filt <- wcvp_accepted_merger %>% 
  filter(is.na(state))


checks <- wcvp_accepted_merger %>% 
  filter(apd_plant_growth_form == "check")


dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% 
  filter(!area == "") %>% 
  filter(!location_doubtful == 1) %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)


data <- wcvp_accepted_merger %>% 
  dplyr::select(lifeform_description, humphreys_lifeform, 
                wood1, state, apd_plant_growth_form, family, genus, plant_name_id) %>% 
  right_join(dist_native, by = "plant_name_id", multiple = "all")

data_ref <- data %>% 
  filter(state == "check") %>% 
  mutate(state = case_when(state == "check" ~ NA)) 



data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) %>% 
  mutate(humphreys_lifeform = ifelse(is.na(humphreys_lifeform), "", humphreys_lifeform), 
         wood1 = ifelse(is.na(wood1), "variable", wood1))




data_test <-  data_merger %>% 
  #filter(family == "Fabaceae") %>% 
  dplyr::select(plant_name_id, state, genus, family, region) #%>% #humphreys_lifeform, wood1,
  #slice_sample(n = 500000)



states <- prodNA(as.data.frame(data_nona$state), noNA = 0.33)
#apd_plant_growth_form <- prodNA(as.data.frame(data_nona$apd_plant_growth_form), noNA = 0.33)

data_nona$state <- states
#data_nona$apd_plant_growth_form <- apd_plant_growth_form

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
data_train_ref <- as.data.frame(unclass(data_nona),stringsAsFactors=TRUE)


mice_rf_imputed <- miceRanger(data_train_ref, 
           vars = c("state"), 
           returnModels = T, 
           m = 5, 
           maxiter = 5, 
           valueSelector = "value", 
           num.trees = 100,
           )


plotDistributions(mice_rf_imputed)
plotCorrelations(mice_rf_imputed)
plotVarConvergence(mice_rf_imputed, vars = "allCategorical")
plotModelError(mice_rf_imputed)
plotVarImportance(mice_rf_imputed, display = "Relative")
plotImputationVariance(mice_rf_imputed)




importance_frame <- data.frame()
for (i in 1: length(mice_rf_imputed$allImps)) {
  data_list <- mice_rf_imputed$allImport[[i]]
  data_frame <- do.call(rbind, data_list)
  data_frame$id <- paste0("dataset_",i)
  importance_frame <- rbind(importance_frame, data_frame)
  
}

importance_frame_rel <- importance_frame %>% 
  mutate(across(where(is.numeric), relative))


error_frame <- data.frame()
for (i in 1: length(mice_rf_imputed$allError)) {
  data_frame <- mice_rf_imputed$allError[[i]]
  data_frame$id <- paste0(i)
  error_frame <- rbind(error_frame, data_frame)
  
}


imps_frame <- data.frame()
for (i in 1: length(mice_rf_imputed$finalImps)) {
  data_list <- mice_rf_imputed$finalImps[i]
  data_frame <- rbindlist(data_list)
  data_frame$id <- paste0(i)
  imps_frame <- rbind(imps_frame, data_frame)
}


imps_frame_final <- data.frame()
for (i in 1: length(mice_rf_imputed$allImps)) {
  data_list <- mice_rf_imputed$allImps[[i]]
  data_frame <- rbindlist(data_list, idcol= T) %>% 
    rename(iteration = .id)
  data_frame$id <- paste0(i)
  imps_frame_final <- rbind(imps_frame_final, data_frame)
}




best <- error_frame %>% 
  slice_max(state, n =1)

imps_multiple <- imps_frame_final %>% 
  group_by(id, iteration,state) %>% 
  reframe(count = length(state))


final_importance <- mice_rf_imputed$finalImport
final_importance_rel <- do.call(rbind, final_importance) %>% 
  mutate(across(where(is.numeric), relative))

#mice_rf_imputed <- paste0("mice_rf_imputed$allImps$Dataset_",best$id,"$Iteration_",best$iteration)
#import <- rbindlist(mice_rf_imputed$allError)
#imput <- as.data.frame(mice_rf_imputed$finalImps$Dataset_1) # highest accuracy
#names(imput) <- "state_imp"


#table(imput$state_imp)

# imputed data
dataList <- completeData(mice_rf_imputed)
data_imputed <- rbindlist(dataList)

length(unique(data_imputed$plant_name_id))
length(unique(data_train_ref$plant_name_id))

# checking for varying states
data_state <- data_imputed %>% 
  group_by(plant_name_id) %>% 
  mutate(count = length(unique(state))) %>%  
  mutate(total = length(plant_name_id)) 

# remove geographic duplicates 
data_state_small <- data_state %>% 
  ungroup() %>%  
  dplyr::select(plant_name_id,  state, count) %>% # lifeform_description
  filter(!duplicated(.))

# filter multistate entries
data_multistate <- data_state %>% 
  filter(count == 2)

# filter unistate entries and remove geographic duplicates
data_unistate_small <- data_state_small %>% 
  filter(count == 1) 

length(unique(data_unistate_small$plant_name_id))

# look into multistate data to determine which ones to keep
data_multistate_exam <- data_multistate %>%
  group_by(plant_name_id, state) %>% 
  reframe(freq = length(state), 
          rel = freq/total) %>% 
  filter(!duplicated(.))

# look at the once with more than 50% in one category
data_multistate_rel <- data_multistate_exam %>% 
  group_by(plant_name_id) %>% 
  slice_max(rel, n = 1) %>% 
  filter(rel > 0.5)


# ectract all 
data_multistate_exam_split <- data_multistate_exam %>% 
  filter(!rel== 0.5)


# association analysis of the two variables
data_na_train <- data_train_ref %>% 
  filter(is.na(state))

data_origin_stat <- data_test %>% 
  filter(plant_name_id %in% data_na_train$plant_name_id) 

data_state_small_real <- data_origin_stat %>% 
  dplyr::select(plant_name_id, state) %>% 
  filter(!duplicated(plant_name_id))

cramerxx <- data_imputed %>% 
  ungroup() %>%  
  dplyr::select(plant_name_id, region, imp_state = state) %>% 
  filter(!duplicated(.)) %>% 
  right_join(data_state_small_real, by = "plant_name_id")


tab <- xtabs(~state + imp_state , data = cramerxx)
summary(assocstats(tab))
mosaicplot(tab)

library(mltools)

actuals <- as.factor(cramerxx$state)
preds <- as.factor(cramerxx$imp_state)
mcc(preds, actuals) # matthwees correlation coefficient


 # mutate(quality = case_when(count == 1 ~ "good",
                 #  count == 2 ~ "check")) 
install.packages("vcd")
data("Arthritis")
tab <- xtabs(~Improved + Treatment, data = Arthritis)
summary(assocstats(tab))
assocstats(UCBAdmissions)

check <- mice_rf_imputed$data
aa <- mice_rf_imputed$naWhere

table(data_test$state)


# model testing



data <- quality %>% 
  dplyr::select(lifeform_description,  humphreys_lifeform, 
                wood1, state, apd_plant_growth_form, quality_check, family, genus, plant_name_id) %>% 
  right_join(dist_native, by = "plant_name_id", multiple = "all")



data_no_na <- na.omit(data)
data_no_na <- as.data.frame(unclass(data_no_na),stringsAsFactors=TRUE) %>% 
  mutate(ID = 1:nrow(data_no_na))

data_test <- data_no_na %>% 
  filter(!quality_check == "check") %>% 
  filter(!state == "check") %>% 
  filter(!apd_plant_growth_form == "check")



data_test2 <- data_test

rows <- sample(1:nrow(data_test), 10000, replace=F)

states <- prodNA(as.data.frame(data_test$state), noNA = 0.33)


data_test$state  <- states$`data_test$state`

data_test_no_na <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
length(unique(data_no_na$ID))

data_test_filt <- data_test_no_na %>% 
  filter(is.na(data_test$state))

data_test_filled <- data_test2 %>% 
  filter(ID %in% data_test_filt$ID)

data_test_ref <-  data_test_no_na %>% 
  #filter(family == "Fabaceae") %>% 
  dplyr::select(ID,  apd_plant_growth_form, state, region, wood1, humphreys_lifeform, family, genus, lifeform_description)

data_test_good <- data_test2 %>% 
  #filter(family == "Fabaceae") %>% 
  dplyr::select(ID,  apd_plant_growth_form, state, region, wood1, humphreys_lifeform, family, genus,lifeform_description)



#imputed_data <- missForest(data_test_ref, xtrue = data_test_good,  maxiter = 10, ntree = 100,  verbose = T)

imputed_data <- missRanger(data_test_ref, formula = state ~ . -ID,  maxiter = 5, num.trees = 500,
                           returnOOB = F, pmm.k	= 0, mtry = 7, verbose = 2, write.fores = T)


imputed_data <- missRanger(data_test_ref, formula = state ~ . -ID,  maxiter = 5, num.trees = 500,
                           returnOOB = F, pmm.k	= 0, mtry = 7, verbose = 2, write.fores = T)


data_test_check <- imputed_data %>%
  filter(ID %in% data_test_filt$ID) %>% 
  dplyr::select(ID, imp_state = state)

data_compare <- data_test_filled %>% 
  filter(ID %in% data_test_check$ID) %>% 
  left_join(data_test_check, by = "ID") %>% 
  dplyr::select(state, imp_state, ID, apd_plant_growth_form, region, wood1, humphreys_lifeform, family, genus)

data_stats <- data_compare %>% 
  mutate(quality = case_when(state == imp_state ~ "good", 
                             !state == imp_state ~ "check"))

data_sumstats <- as.data.frame(table(data_stats$quality))
data_sumstats$Prop <- data_sumstats$Freq / sum(data_sumstats$Freq)
data_sumstats

