library(DEoptim)
library(progress)
library(tidyverse)
library(data.table)
library(sf)
library(foreach)
library(doSNOW)



output <- function(species=NULL,best_corr=NULL) {
  me <- data.frame(
    Species = species,
    best_corr = best_corr
  )
  
  return(me)
}

objective <- function(subset_indices) {
  species <- slice_sample(n = round(subset_indices),plantlist_names)
  dist <- plantlist_dist %>%
    filter(plant_name_id %in% species$plant_name_id)
  sample_rich_bru <- dist %>%
    group_by(area_code_l3) %>%
    summarise(richness_sample = n_distinct(plant_name_id)) %>%
    rename(LEVEL_COD = area_code_l3)
  rich_rel <- rich_overall_bru %>%
    left_join(sample_rich_bru, by = "LEVEL_COD") %>%
    replace(is.na(.), 0) %>%
    mutate(sp = length(species))
  cor.sp <- cor.test(rich_rel$richness_sample, rich_rel$richness, 
                     method = "spearman", exact = FALSE)[[4]]
  return(-cor.sp)  # Return negative correlation (since DEoptim minimizes)
  #return(ifelse(correlation >= 0.95, correlation, 0))
}


# Working directory ------------------------------------------------------------
# Import data ------------------------------------------------------------------
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

richness_patterns <- fread("data/richness_patterns.txt")
plants_full <- fread("data/wcvp_accepted_merged.txt")
dist_native <- fread("data/dist_native.txt")


sample <- sample(plants_full$plant_name_id, 10000)

plants_sample <- plants_full %>%  
  filter(plant_name_id %in% sample)

dist_sample <- dist_native %>%  
  filter(plant_name_id %in% sample)

# creating objects for the function

plantlist_names <- plants_full %>% 
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)




richness_patterns_con <- plantlist_dist %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- plantlist_dist %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- plantlist_dist %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))








#clusterExport(cl = NULL, varlist, envir = .GlobalEnv)
# Ver 1

# Define lower and upper bounds for each species (dimension)
lower_bound <- 1  # Start from the first species
upper_bound <- nrow(plantlist_names) # Go up to the last species

# Initialize DE parameters


results_list <- vector("list", length = 1390)

# Define the number of cores you want to use
num_cores <- 60  # Adjust this based on your system

# Create a cluster
cl <- makeCluster(num_cores)


# Register the cluster with foreach
registerDoSNOW(cl)



iterations <- 5123
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


solutions <- foreach(i = 1:5123, 
                     .packages = c("tidyverse", "DEoptim"), 
                     .combine = "rbind",  
                     .options.snow = opts,
                     .verbose = T) %dopar% 
  {
                       
    de_result <- DEoptim(objective, lower = i, upper = i,
                         DEoptim.control(VTR = -Inf, strategy = 2,
                                         bs = FALSE, itermax = 10000, CR = 0.5, F = 0.8, 
                                         storepopfreq = 1, trace = FALSE, parallelType = "foreach"))
    
    output(de_result$optim$bestmem, de_result$optim$bestval*-1)
    }

write.table(solutions, "data/optimised/DEoptim_results.txt")
close(pb)
stopCluster(cl)


