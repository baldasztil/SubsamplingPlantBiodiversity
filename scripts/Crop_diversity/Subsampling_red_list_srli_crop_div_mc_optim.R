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


output <- function(species=NULL,best_corr=NULL)
{
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


full_redlist <- read.csv("data/red/redlist_full_syn.csv")
redlist_index <- read.csv("data/red/srli_full.csv")
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
dist_native <- dist_native_all %>% 
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


results_list <- vector("list", length = 1390)

# Define the number of cores you want to use
num_cores <- 30  # Adjust this based on your system

# Create a cluster
cl <- makeCluster(num_cores)


# Register the cluster with foreach
registerDoSNOW(cl)



solutions <- foreach(i = 1:5123, 
                     .packages = c("tidyverse", "DEoptim"), 
                     .combine = "rbind",  
                     .options.snow = opts, 
                     .verbose = T) %dopar% 
  {
    
    de_result <- DEoptim(objective, lower = i, upper = i,
                         DEoptim.control(VTR = -Inf, strategy = 2,
                                         bs = FALSE, itermax = 1000, CR = 0.5, F = 0.8, 
                                         storepopfreq = 1, trace = FALSE, parallelType = "foreach"))
    write.table(i, file=paste("data/optimised/outfile/outfile_",i,"_Species_corr",round(de_result$optim$bestval*-1, digits = 3),".txt"))
    output(de_result$optim$bestmem, de_result$optim$bestval*-1)
  }

write.table(solutions, "data/optimised/DEoptim_results_red_list.txt")
stopCluster(cl)

# SRLI -------------------------------------------------------------------------

plants_full <- redlist_index %>% 
  filter(plant_name_id %in% redlist_index$plant_name_id) 


dist_native <- dist_native_all %>% 
  filter(plant_name_id %in% redlist_index$plant_name_id) 

#### looping 
# names 
plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

# distribution
plantlist_dist <- dist_native %>% 
  dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

# Region names 
tdwg_codes <- as.data.frame(tdwg_3) %>% 
  dplyr::select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

# Overall patterns at botanical recording unit scale  
rich_overall_bru <- richness_patterns_allplants %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))


#### Analysis 
for (i in 1:100)  {
  samples <- list()
  
  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(1,nrow(plantlist_names),1), mc.cores = 30, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/red/srli/Samples_iteration_srli_",i,".txt"))
}
