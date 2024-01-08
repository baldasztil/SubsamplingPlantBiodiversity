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
library(parallel)
library(doParallel)
library(sf)
options(dplyr.summarise.inform = FALSE)
# Defining functions -----------------------------------------------------------


# this is a function to create a subsample the WCVP and calculate correlations
# to the global patterns of taxonomic richness at a family level
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
  # this creates a subsample
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
    filter(plantlist_names$plant_name_id %in% cumulative_namelist)
  
  dist <- dist_native %>% filter(plant_name_id %in% cumulative_namelist)
  
  # richness patterns across TDWG Level 3 region split
  #for each family within a region
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3, family) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  # combining with the overall patterns of taxonomic richness
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = c("LEVEL_COD", "family")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist)) %>% 
    ungroup()
  
  rich_ex <- sample_rich_bru %>%  
    filter(!family %in% rich_rel$family)
  
  cor.test(rich_rel$richness_sample, rich_rel$richness,
                                method="spearman", exact =F)
  
  # measuring correlation to the overall patterns of taxonomic richness
  rich_rel_bru <- rich_rel %>%
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent 
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD, .family) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD, .family) %>% 
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
  
  # extracting cumulative patterns  
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}

# this is a for loop around the function
parallel.subsampling <- function (iteration, sample_n) {
  for (i in 1:iteration)  {
    samples <- list()
    
    print(paste0("This is iteration ", i)) 
    samples[[1]] <- mclapply(seq(0,nrow(plantlist_names),sample_n), mc.cores = 30, subsampling.plants)
    xx <- do.call(bind_rows, samples[[1]])
    write.table(xx, paste0("data/taxonomic_rep_5more/Samples_iteration_",i,".txt"))
  }
  
  return(samples)
}

# Import data ------------------------------------------------------------------

tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

plants_full <- fread("data/wcvp_accepted_merged.txt")

plant_family <- plants_full %>% 
  select(plant_name_id, family)
dist <- fread("data/dist_native.txt")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)


# baseline richness ------------------------------------------------------------

# calculating overall taxonomic richness patterns at family level

dist_native <- dist %>% 
  left_join(plant_family, by = "plant_name_id")


richness_patterns_con <- dist_native %>% 
  group_by(continent_code_l1, family) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_con$LEVEL1_COD <- as.character(richness_patterns_con$LEVEL1_COD)

richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2,family) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")
richness_patterns_reg$LEVEL2_COD <- as.character(richness_patterns_reg$LEVEL2_COD)

richness_patterns_bru <- dist_native %>% 
  group_by(area_code_l3,family) %>% 
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
                           richness_patterns3) %>% 
  ungroup() %>% 
  mutate(ID2 = row_number())

length(unique(richness_patterns$family))

# Analysis ---------------------------------------------------------------------

# creating objects for function
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD")) %>% 
  filter(richness > 5) #%>% 
  
  #group_by(LEVEL_COD) %>% 
  #slice_max(richness, n = 20)
rich_ex <- richness_patterns %>% 
  filter(!ID2 %in% rich_overall_bru$ID2)

###### running the function

parallel.subsampling(100, 1)

