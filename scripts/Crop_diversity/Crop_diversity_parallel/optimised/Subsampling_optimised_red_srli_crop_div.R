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
library(sf)
library(parallel)
library(DEoptim)
library(foreach)
library(doSNOW)
library(progress)


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


full_redlist <- read.csv("data/red/srli_full.csv")
richness_patterns_allplants <- fread("data/richness_patterns.txt")



plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")


tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)



# Subsampling ------------------------------------------------------------------

# Redlist ----------------------------------------------------------------------

# names 
plantlist_names <- full_redlist %>%  
   dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

# distribution
plantlist_dist <- dist_native_all %>% 
  filter(plant_name_id %in% full_redlist$plant_name_id)

# Region names 
tdwg_codes <- as.data.frame(tdwg_3) %>% 
   dplyr::select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

# Overall patterns at botanical recording unit scale  
rich_overall_bru <- richness_patterns_allplants %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))


# Define lower and upper bounds for each species (dimension)
lower_bound <- 1  # Start from the first species
upper_bound <- nrow(plantlist_names) # Go up to the last species

# Initialize DE parameters

# Define the number of cores you want to use
num_cores <- 30  # Adjust this based on your system

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
                                         bs = FALSE, itermax = 100, CR = 0.5, F = 0.8, 
                                         storepopfreq = 1, trace = FALSE, parallelType = "foreach"))
    output(de_result$optim$bestmem, de_result$optim$bestval*-1)
  }

write.table(solutions, "data/optimised/DEoptim_results_srli_100.txt")
stopCluster(cl)
