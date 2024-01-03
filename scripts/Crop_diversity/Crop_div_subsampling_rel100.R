
# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------




# Libraries --------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(purrr)
library(progress)
library(data.table)
library(ggplot2)
# Defining functions -----------------------------------------------------------

# Working directory ------------------------------------------------------------
# Import data ------------------------------------------------------------------

wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
tdwg_codes <- fread("data/tdwg_codes.csv", header = T) 

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
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)

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
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")


richness_patterns_bru <- dist_native %>% 
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

#-------------------------------------------------------------------------------


# Adjust data ------------------------------------------------------------------

#table of number of unique species per botanical unit here at different levels 


# Analysis ---------------------------------------------------------------------
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

#plantlist_names <- plants_full %>%  
 # select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name) %>%
  #sample_n(350250)

rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

######
start.time <- Sys.time()
cumulative_namelist <- data.frame()
rich_rel_cumulative <- data.frame()
rich_cumulative <- data.frame()
samplelist <- list()
samples <- list()
for (x in 1:100) {
    repeat{
      # randomly sampling the data for species ids 
      plantlist_names_left <- plantlist_names %>% 
        filter(!plant_name_id %in%  cumulative_namelist$plant_name_id)
      # too include last sample (smaller than 500)
      if (nrow(plantlist_names_left) > 100) {
        species_sample <- sample_n(plantlist_names_left, 100)
      } else {
        species_sample <- plantlist_names_left
      }
      cumulative_namelist <- rbind(cumulative_namelist, species_sample)
      species <- cumulative_namelist
      
      # caluclating sample species richness
      dist <- dist_native %>% 
        filter(plant_name_id %in% species$plant_name_id)
      
      # richness patterns across brus
      sample_rich_bru <- dist %>% 
        group_by(area_code_l3) %>% 
        summarise(richness_sample = n_distinct(plant_name_id)) %>% 
        rename(LEVEL_COD = area_code_l3)
      
      rich_rel <- rich_overall_bru %>% 
        left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
        replace(is.na(.), 0) %>% 
        mutate(sp = nrow(cumulative_namelist))
      # measuring correlation
      
      # overall
      rich_rel_bru <- rich_rel %>% 
        summarise(cor.sp = cor.test(richness_sample, richness,
                                               method="spearman", exact =F)[[4]],
                  cor.pe = as.numeric(cor(richness_sample, richness, 
                                               method="pearson"))) %>% 
        mutate(id = "overall",
               sp = nrow(cumulative_namelist))
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
               sp = nrow(cumulative_namelist)) %>% 
        left_join(spear, by = "id")
  
        # cumulative patterns
        rich_rel_cumulative <- rbind(rich_rel_cumulative,rich_rel_bru, rich_rel_con)
      if (nrow(cumulative_namelist) == nrow(plantlist_names)){
        break
      }
    if (nrow(cumulative_namelist) %in% seq(0,nrow(plantlist_names),1000)){
    print(paste0("There are ", nrow(cumulative_namelist) ," species in the subsample")) 
    }
  }
    samplelist[[x]] <- rich_rel_cumulative 
    cumulative_namelist <- data.frame()
    rich_rel_cumulative <- data.frame()
    print(paste0("This is iteration ", x))
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
}


data_out <- as.data.frame(do.call(rbind, samplelist))
write.table(data_out, "full_samples_rel.txt")
