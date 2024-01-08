# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(parallel)
# Defining functions -----------------------------------------------------------
subsampling.plants <- function(plant_family) {
  plantlist_names <- plantlist %>% 
    filter(family %in% plant_family)
  list_rich_rel_cumulative <- list()
  for (i in seq(1,nrow(plantlist_names),1)) {
    cumulative_namelist <- c()
    plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist)
  if (nrow(plantlist_names_left) > i) {
    species_sample <- sample_n(plantlist_names_left, i)
    # too include last sample
  } else {
    species_sample <- plantlist_names_left
  }
  cumulative_namelist <- c(cumulative_namelist, species_sample$plant_name_id)
  species <-plantlist_names %>% 
    filter(plant_name_id %in%  cumulative_namelist)
  species_names <- species$plant_name_id
  dist <- dist_native %>% 
    filter(plant_name_id %in% species_names)
  
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  rich_rel$family <- unique(species$family)
  # measuring correlation
  # overall
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(richness_sample, richness,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(richness_sample, richness,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  rich_rel_bru$family <- unique(species$family)
  # per continent 
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$richness_sample, .$richness, 
                      method="pearson", data = ., exact =F)[[4]])  
  # bind all
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  rich_rel_con$family <- unique(species$family)
  
  list_rich_rel_cumulative[[i]] <- bind_rows(rich_rel_bru, rich_rel_con)
  #if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
  # print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  #}
  
  # cumulative pattern
  }
  print(paste0("This is ", plant_family)) 
  return(as.data.frame(do.call(bind_rows, list_rich_rel_cumulative)))
}
parallel.subsampling <- function (iteration, plant_families) {
  
  for (i in 1:iteration)  {
    samples <- list()
    
    print(paste0("This is iteration ", i))
    samples[[1]] <- mclapply(plant_families, mc.cores = 30, subsampling.plants)
    xx <- as.data.frame(do.call(bind_rows, samples))
    write.table(xx, paste0("data/families_rerun/Samples_iteration_families_",i,".txt"))
  }
   
  return(samples)
}



#family_list <- plantlist %>%
 # mutate(family2 = family) %>% 
  #group_by(family2) %>% 
  #group_map(~.x)


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

plant_families <- unique(plants_full$family)
# analyis ----------------------------------------------------------------------

plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))


samples <- parallel.subsampling(100, plant_families)



family_list_check <- data_out %>%
  group_by(family) %>%
  summarise(species = max(sp)) %>%
  left_join(family_list, by ="family") %>% 
  mutate(diff = species.x - species.y)

unique(family_list_check$diff)