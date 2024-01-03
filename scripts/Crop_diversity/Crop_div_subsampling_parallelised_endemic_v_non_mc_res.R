# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------

# This is a script to create random subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species richness in the subsample to global patterns. We split the 
# dataset into endemic and non endemic species


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(ggplot2)
library(parallel)
# Defining functions -----------------------------------------------------------

# this is a function to subsample the WCVP
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
  # accumulative subsampling 
  plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist)
  if (nrow(plantlist_names_left) > spec_n) {
    species_sample <- sample_n(plantlist_names_left, spec_n)
    # too include last sample
  } else {
    species_sample <- plantlist_names_left
  }
  cumulative_namelist <- c(cumulative_namelist, species_sample$plant_name_id)
  species <-plantlist_names %>% 
    filter(plant_name_id %in%  cumulative_namelist)
  
  dist <- dist_native %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # richness patterns across tdwg3 areas
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  # measuring correlation
  
  #correlation to overall richness 
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent correlation
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="pearson", data = ., exact =F)[[4]])  
  
  # combining the results
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  
  if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
    print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  
  #extracting cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
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


# manipulate data --------------------------------------------------------------
# extracting names from shapefile  
tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

# making it a character for later
tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

tdwg_codes <- as.data.frame(tdwg_3) %>% 
  select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

tdwg_codes$LEVEL1_COD <- as.character(tdwg_codes$LEVEL1_COD)

index_endemic <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n == 1) # 56%

index_widespread <- dist_native %>% 
  group_by(plant_name_id) %>% 
  summarise(n = n()) %>% 
  filter(n > 1) # 43.9%

# endemic ----------------------------------------------------------------------

# creating objects for the function
plantlist_dist <- dist_native %>% 
  filter(plant_name_id %in% index_endemic$plant_name_id) %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>% 
  filter(plant_name_id %in% index_endemic$plant_name_id) %>% 
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))


#### Analysis 
for (i in 1:100)  {
  samples <- list()

  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(0,nrow(plantlist_names),1), mc.cores = 30, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/end_nonend/Samples_iteration_end_",i,".txt"))
}


# non-endemic ------------------------------------------------------------------

# creating objects for the function
plantlist_dist <- dist_native %>% 
  filter(plant_name_id %in% index_widespread$plant_name_id) %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>% 
  filter(plant_name_id %in% index_widespread$plant_name_id) %>% 
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

#### Analysis 
for (i in 1:100)  {
  samples <- list()
  
  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(0,nrow(plantlist_names),1), mc.cores = 30, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/end_nonend/Samples_iteration_nonend_",i,".txt"))
}


