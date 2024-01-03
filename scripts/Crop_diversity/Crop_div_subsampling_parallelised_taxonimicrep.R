# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(parallel)
library(doParallel)
library(sf)
# Defining functions -----------------------------------------------------------
subsampling.plants <- function(spec_n) {
  cumulative_namelist <- c()
  list_rich_rel_cumulative <- list()
  
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
  
  # richness patterns across brus
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3, family) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3)
  
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = c("LEVEL_COD", "family")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist)) %>% 
    ungroup()
  # measuring correlation
  
  # overall
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
  
  # bind all
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  
  #if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
  # print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  #}
  
  # cumulative pattern
  if (length(cumulative_namelist) %in% seq(0,nrow(plantlist_names),10000)){
  print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  return(bind_rows(rich_rel_bru, rich_rel_con))
}
parallel.subsampling <- function (iteration,sample_n) {
  samples <- list()
  for (i in 1:iteration)  {
    print(paste0("This is iteration ", iteration))
    samples[[i]] <- parLapply(clust, seq(0,nrow(plantlist_names),sample_n), subsampling.plants)
  }
  return(samples)
}

# Working directory ------------------------------------------------------------

# Import data ------------------------------------------------------------------

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

# manipulate data --------------------------------------------------------------

wcvp_placed <- wcvp_raw %>% 
  filter(!accepted_plant_name_id %in% c(""))

wcvp_sp <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(!taxon_rank %in% c("Genus"))

wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) 


dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% 
  filter(!area == "") %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id) %>% 
  left_join(wcvp_accepted, by = "plant_name_id")

dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))


plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")

plants_full_extinct_norange <- wcvp_accepted %>% 
  filter(!plant_name_id %in% plants_full$plant_name_id)

# baseline richness ------------------------------------------------------------
# calculating overall richness patterns 

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
                           richness_patterns3)
#-------------------------------------------------------------------------------


# Adjust data ------------------------------------------------------------------

#table of number of unique species per botanical unit here at different levels 


# Analysis ---------------------------------------------------------------------
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

######
clust <- makeCluster(6, outfile = "log_taxrich.txt")
clusterExport(clust, c("dist_native","plantlist_dist","plantlist_names","rich_overall_bru"))
clusterEvalQ(clust, library(tidyverse))
doParallel::registerDoParallel(clust)

samples_rel <- parallel.subsampling(100, 1)
stopCluster(clust)

data_out <- as.data.frame(do.call(bind_rows, samples_rel))
write.table(data_out, "families_rel_stats.txt")
