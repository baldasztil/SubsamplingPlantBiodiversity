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
summary.stats <- function(x, v1) {
  true_error_frame <- data.frame()
  for (i in (1:x)) {
    
    v1$plant_name_id <- data_test_na$plant_name_id
    v1$apd_plant_growth_form_origin <- data_test_na$apd_plant_growth_form
    
    a <- v1 %>% 
      dplyr::select(ends_with(paste0(i)), plant_name_id, apd_plant_growth_form_origin)
    names(a) <- c("apd_plant_growth_form_imp", "plant_name_id", "apd_plant_growth_form_origin")
    
    b <-   as.data.frame(xtabs(~ apd_plant_growth_form_imp + apd_plant_growth_form_origin, data = a)) %>% 
      mutate(apd_plant_growth_form_origin = as.character(apd_plant_growth_form_origin), 
             apd_plant_growth_form_imp = as.character(apd_plant_growth_form_imp)) %>% 
      mutate(true = apd_plant_growth_form_imp == apd_plant_growth_form_origin)  %>% 
      group_by(apd_plant_growth_form_origin) %>% 
      mutate(total = sum (Freq))
    
    c <- b %>%
      filter(true == FALSE) %>%
      mutate(rel = round(Freq / total , 3))
    #plot(b$rel ~ b$apd_plant_growth_form_origin)
    
    d <-  c %>% 
      ungroup() %>% 
      group_by(apd_plant_growth_form_origin) %>% 
      reframe(rel_error = round(sum(Freq)/ total, 3)) %>% 
      filter(!duplicated(.))
    
    e <- b %>% 
      ungroup() %>% 
      group_by(true) %>% 
      reframe(sum = sum(Freq)) %>% 
      ungroup() %>% 
      mutate(total = sum(sum),
             rel_error = round(sum/total, 3))
    
    f <- e %>% 
      dplyr::select(apd_plant_growth_form_origin = true, rel_error) %>% 
      mutate(apd_plant_growth_form_origin = as.character(apd_plant_growth_form_origin))
    
    g <- rbind(d, f)
    g$dataset <- i 
    
    true_error_frame <- rbind(g, true_error_frame)
  }
  true_error_frame_means <- true_error_frame %>% 
    group_by(apd_plant_growth_form_origin) %>% 
    reframe(mean_error = mean(rel_error), 
            sd = sd(rel_error))
  return(true_error_frame_means)
  
}

# Import data ------------------------------------------------------------------
plants_full <- fread("data/wcvp_accepted_merged.txt") 
dist_native <- fread("data/dist_native.txt") 
#woodiness <- read.table("data/traits_all/woodiness_wcvp_data.txt", header = T, fill = T)  

#woodiness <- read.csv("data/traits_all/woodi_herbaceous_rf.csv")

lifeform <- read.csv("data/traits_all/lifeform_mapping_original.csv") %>% 
  dplyr::select(lifeform_description, apd_plant_growth_form) %>% 
  mutate(apd_plant_growth_form = ifelse(apd_plant_growth_form =="subshrub", "shrub",
                                        apd_plant_growth_form)) %>% 
  mutate(apd_plant_growth_form = gsub(" ", "", apd_plant_growth_form))

plants_old_growth <- fread("data/wcvp/wcvp_names.txt") %>% 
  left_join(lifeform, by = "lifeform_description") %>% 
  filter(taxon_name %in% plants_full$taxon_name) %>% 
  dplyr::select(taxon_name, apd_plant_growth_form, old_lifeform = lifeform_description) %>% 
  filter(!duplicated(taxon_name))



# manipulate data --------------------------------------------------------------

wcvp_accepted_merger <- plants_full %>%   
  left_join(plants_old_growth, by = "taxon_name") 

growth_summary <- wcvp_accepted_merger %>% 
  group_by(lifeform_description, apd_plant_growth_form, old_lifeform) %>% 
  summarise() 

growth_summary2 <- wcvp_accepted_merger %>% 
  group_by(lifeform_description, apd_plant_growth_form) %>% 
  summarise() 


# manipulate data --------------------------------------------------------------


data <- wcvp_accepted_merger %>% 
  dplyr::select(lifeform_description, 
                apd_plant_growth_form, family, genus, plant_name_id, continent)

data_ref <- data %>% 
  filter(apd_plant_growth_form == "check" | is.na(apd_plant_growth_form)) %>% 
  mutate(apd_plant_growth_form = NA)

data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) 

data_test <-  data_merger %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v1 <- miceRanger(data_test_ref, 
                                 vars = c("apd_plant_growth_form"), 
                                 returnModels = T, 
                                 m = 99, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 501
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
  data_frame <- as.data.frame(data_list$Dataset$apd_plant_growth_form) %>% 
    rename(state = "data_list$Dataset$apd_plant_growth_form")
  data_frame$id <- paste0(i)
  imps_frame <- rbind(imps_frame, data_frame)
}



best <- error_frame %>% 
  slice_max(apd_plant_growth_form, n =1)


final_importance <- mice_rf_imputed_v1$finalImport
final_importance_rel <- do.call(rbind, final_importance) %>% 
  mutate_all(funs(ifelse(is.na(.), 0, .))) %>%  
  mutate(across(where(is.numeric), relative))

# imputed data
dataList <- completeData(mice_rf_imputed_v1)
data_imputed <- rbindlist(dataList)
length(unique(data_imputed$plant_name_id))

nas <- data_test_ref %>% 
  filter(is.na(apd_plant_growth_form))

imputed_apd_plant_growth_formsv1 <- data_imputed %>% 
  filter(plant_name_id %in% nas$plant_name_id) %>% 
  dplyr::select(plant_name_id_v1 = plant_name_id, apd_plant_growth_forms_v1 = apd_plant_growth_form)


length(unique(data_imputed$plant_name_id))
length(unique(data_test_ref$plant_name_id))

# checking for varying apd_plant_growth_forms
data_apd_plant_growth_form <- data_imputed %>% 
  group_by(plant_name_id) %>% 
  mutate(count = length(unique(apd_plant_growth_form)), 
         total = length(plant_name_id))

# filter multiapd_plant_growth_form entries
data_multiapd_plant_growth_form <- data_apd_plant_growth_form %>% 
  filter(count >= 2)

length(unique(data_multiapd_plant_growth_form$plant_name_id))

# filter uniapd_plant_growth_form entries and remove geographic duplicates
data_uniapd_plant_growth_form_small <- data_apd_plant_growth_form %>% 
  filter(count == 1) 

length(unique(data_uniapd_plant_growth_form_small$plant_name_id))

# look into multiapd_plant_growth_form data to determine which ones to keep
data_multiapd_plant_growth_form_exam <- data_apd_plant_growth_form %>%
  ungroup() %>% 
  group_by(plant_name_id, apd_plant_growth_form) %>% 
  mutate(freq = length(plant_name_id), 
         rel = freq / total) %>% 
  ungroup() %>% 
  filter(!duplicated(.))


length(unique(data_multiapd_plant_growth_form_exam$plant_name_id))

# look at the once with more than 50% in one category
data_multiapd_plant_growth_form_rel_adj <- data_multiapd_plant_growth_form_exam %>% 
  group_by(plant_name_id) %>% 
  slice_max(rel, n = 1, with_ties = F) %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, certainty = rel) %>% 
  ungroup()

length(unique(data_multiapd_plant_growth_form_rel_adj$plant_name_id))


aa <- data_multiapd_plant_growth_form_exam %>% 
  filter(!plant_name_id %in% data_multiapd_plant_growth_form_rel_adj$plant_name_id)

 

save.image("RF_99_datasets_growths_apd.RData")
write.csv(data_multiapd_plant_growth_form_rel_adj, "data/traits_all/growths_rf.csv")



