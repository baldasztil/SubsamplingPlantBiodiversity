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
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
#woodiness <- read.table("data/traits_all/woodiness_wcvp_data.txt", header = T, fill = T)  

woodiness <- read.csv("data/traits_all/woodi_herbaceous_rf.csv")

lifeform <- read.csv("data/traits_all/lifeform_mapping_original.csv") %>% 
  dplyr::select(lifeform_description, apd_plant_growth_form)

growthform <- fread("data/traits_all/traits_all.csv")
table(growthform$growth_form)
table(lifeform$apd_plant_growth_form)

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
  left_join(lifeform, by = "lifeform_description") %>% 
  left_join(woodiness, by = "plant_name_id") %>% 
  mutate(lifeform_description = ifelse(lifeform_description == "", NA, lifeform_description))

wcvp_accepted_merger_filt <- wcvp_accepted_merger %>% 
  filter(is.na(apd_plant_growth_form))

growth_gift <- growthform %>% 
  dplyr::select(work_species, growth_form) %>% 
  left_join(wcvp_accepted_merger, by = c("work_species" = "taxon_name")) %>% 
  filter(!duplicated(.))

dups_growth_gift_names <- growth_gift %>% 
  filter(duplicated(plant_name_id))

dups_growth_gift <- growth_gift %>% 
  filter(plant_name_id %in% dups_growth_gift_names$plant_name_id)


growth_gift_var <- growth_gift %>% 
  mutate(growth_form = ifelse(plant_name_id %in% dups_growth_gift$plant_name_id, "check", growth_form)) %>% 
  filter(!duplicated(.) & !is.na(plant_name_id))  

length(unique(growth_gift_var$plant_name_id))

xtabs(~ apd_plant_growth_form + growth_form, data = growth_gift_var)

growth_gift_join <- growth_gift_var%>% 
  dplyr::select(plant_name_id, growth_gift = growth_form)

xtabs(~ apd_plant_growth_form + state, data = wcvp_accepted_merger)


data <- wcvp_accepted_merger %>% 
  dplyr::select(lifeform_description, 
                state, apd_plant_growth_form, family, genus, plant_name_id, continent) %>% 
  left_join(growth_gift_join, by = "plant_name_id")

data_ref <- data %>% 
  filter(apd_plant_growth_form == "check" | apd_plant_growth_form == "check " | is.na(apd_plant_growth_form)) %>% 
  mutate(apd_plant_growth_form = NA, 
         lifeform_description = ifelse(lifeform_description == "",NA, lifeform_description)) 


data_merger <- data %>%  
  filter(!plant_name_id %in% data_ref$plant_name_id) %>% 
  rbind(data_ref) 


check_na_in(data_merger)

aa <- as.data.frame(check_na_in(data_ref, true_total = T))
names(aa) <- "missing"
aa$total <- nrow(data_ref)

aa$rel <- aa$missing / aa$total

checks_life_apd <- data_ref %>% 
  filter(is.na(lifeform_description) & is.na(apd_plant_growth_form))

checks_life <- data_ref %>% 
  filter(is.na(lifeform_description))

checks_growth <- data_ref %>% 
  filter(is.na(growth_gift))

checks_growth_check <- data_ref %>% 
  filter(growth_gift == "check")

table(checks_growth$growth_gift)

check_na_in(data_ref)

data_nona  <- data_merger %>% 
  filter(!is.na(apd_plant_growth_form)) %>% 
  arrange(plant_name_id)


data_na_toinsert <- data_nona %>% 
  slice_sample(prop = nrow(data_ref) / nrow(data_merger)) %>% 
  arrange(plant_name_id)

na_life <- round(nrow(data_na_toinsert) * (nrow(checks_life) / nrow(data_ref)))
na_growth <- round(nrow(data_na_toinsert) * (nrow(checks_growth) / nrow(data_ref)))
na_growth_check <- round(nrow(data_na_toinsert) * (nrow(checks_growth_check) / nrow(data_ref))) * 0.75



data_na <- data_na_toinsert %>%  
  mutate(apd_plant_growth_form = NA)
 
table(data_ref$growth_gift)
table(data_nona$growth_gift)


data_na_life <- data_na %>% 
  slice_sample(n = na_life)

data_na_life_gift <- data_na_life %>% 
  slice_sample(n = na_growth)

data_check_gift <- data_na_toinsert %>% 
  filter(!plant_name_id %in% data_na_life$plant_name_id) %>% 
  slice_sample(n = na_growth_check)

data_na_Life_gift_ref <- data_na %>% 
  mutate(lifeform_description = ifelse(plant_name_id %in% data_na_life$plant_name_id, NA, lifeform_description), 
         growth_gift = ifelse(plant_name_id %in% data_na_life_gift$plant_name_id , NA, growth_gift), 
         growth_gift = ifelse(plant_name_id %in% data_check_gift$plant_name_id , "check", growth_gift))

table(data_na_Life_gift_ref$growth_gift)
check_na_in(data_na_Life_gift_ref)

data_na_full <- data_na %>% 
  filter(!plant_name_id %in% data_na_Life_gift_ref$plant_name_id) %>% 
  bind_rows(data_na_Life_gift_ref) %>% 
  arrange(plant_name_id)


check_na_in(data_na_full)
check_na_in(data_ref)


data_sim <- data_nona %>% 
  filter(!plant_name_id %in% data_na_full$plant_name_id) %>% 
  rbind(data_na_full) %>% 
  arrange(plant_name_id)

check_na_in(data_sim)
check_na_in(data_merger)


data_origin <- data_nona

table(data_sim$apd_plant_growth_form)



data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family, lifeform_description, state, growth_gift)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
str(data_test_ref)

mice_rf_imputed_v1 <- miceRanger(data_test_ref, 
                                 vars = "apd_plant_growth_form", 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 2000,
)

data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family, lifeform_description, state, growth_gift) 
data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
str(data_test_ref)

mice_rf_imputed_v2 <- miceRanger(data_test_ref, 
                                 vars = c("lifeform_description", "growth_gift","apd_plant_growth_form"),  
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 2000,
)




data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family, state)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v3 <- miceRanger(data_test_ref, 
                                 vars = c("apd_plant_growth_form"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 1500,
)

data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family)

data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)

mice_rf_imputed_v4 <- miceRanger(data_test_ref, 
                                 vars = c("apd_plant_growth_form"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 1500,
)



data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family, state, growth_gift)
data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
str(data_test_ref)
check_na_in(data_test)
mice_rf_imputed_v5 <- miceRanger(data_test_ref, 
                                 vars = c("growth_gift", "apd_plant_growth_form"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 2000,
)


data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family,  state, growth_gift) %>% 
  mutate(growth_gift = ifelse(growth_gift == "check", NA, growth_gift))
data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
str(data_test_ref)
check_na_in(data_test)

mice_rf_imputed_v6 <- miceRanger(data_test_ref, 
                                 vars = c("growth_gift", "apd_plant_growth_form"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 2000,
)


data_test <-  data_sim %>% 
  dplyr::select(plant_name_id, apd_plant_growth_form, genus, continent, family,  state, lifeform_description)
data_test_ref <- as.data.frame(unclass(data_test),stringsAsFactors=TRUE)
str(data_test_ref)
check_na_in(data_test)

mice_rf_imputed_v7 <- miceRanger(data_test_ref, 
                                 vars = c("lifeform_description", "apd_plant_growth_form"), 
                                 returnModels = T, 
                                 m = 3, 
                                 maxiter = 5, 
                                 valueSelector = "value", 
                                 num.trees = 2000,
)





nas <- mice_rf_imputed_v1$data %>% 
  filter(is.na(apd_plant_growth_form))

data_test_na <- data_origin %>% 
  filter(plant_name_id %in% nas$plant_name_id) %>% 
  arrange(plant_name_id)



v1 <- as.data.frame(getVarImps(mice_rf_imputed_v1, 1:3, var = "apd_plant_growth_form"))
perf_v1 <- summary.stats(3, v1)
perf_v1$v <- "v1"

v2 <- as.data.frame(getVarImps(mice_rf_imputed_v2, 1:3, var = "apd_plant_growth_form"))
perf_v2 <- summary.stats(3, v2)
perf_v2$v <- "v2"


v3 <- as.data.frame(getVarImps(mice_rf_imputed_v3, 1:3, var = "apd_plant_growth_form"))
perf_v3 <- summary.stats(3, v3)
perf_v3$v <- "v3" 

v4 <- as.data.frame(getVarImps(mice_rf_imputed_v4, 1:3, var = "apd_plant_growth_form"))
perf_v4 <- summary.stats(3, v4)
perf_v4$v <- "v4"  

v5 <- as.data.frame(getVarImps(mice_rf_imputed_v5, 1:3, var = "apd_plant_growth_form"))
perf_v5 <- summary.stats(3, v5)
perf_v5$v <- "v5"

v6 <- as.data.frame(getVarImps(mice_rf_imputed_v6, 1:3, var = "apd_plant_growth_form"))
perf_v6 <- summary.stats(3, v6)
perf_v6$v <- "v6"

#v7 <- as.data.frame(getVarImps(mice_rf_imputed_v7, 1:3, var = "apd_plant_growth_form"))
#perf_v7 <- summary.stats(3, v7)
#perf_v7$v <- "v7"


perf <- rbind(perf_v1, perf_v2, perf_v3, perf_v4, perf_v5, perf_v6)

perf_stats <- perf %>% 
  filter(apd_plant_growth_form_origin %in% c("TRUE")) %>% 
  group_by(apd_plant_growth_form_origin) %>% 
  slice_max(mean_error, n = 5)

perf_v3

plotDistributions(mice_rf_imputed_v3)
plotCorrelations(mice_rf_imputed_v3)
plotVarConvergence(mice_rf_imputed_v3, vars = "allCategorical")
plotModelError(mice_rf_imputed_v3)
plotVarImportance(mice_rf_imputed_v3, display = "Relative")
plotImputationVariance(mice_rf_imputed_v3)


summary(assocstats(tab))
mcc(v1_imp$state_v1, v2_imp$state_v2 ) # matthwees correlation coefficient


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


