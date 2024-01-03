# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------
  # This is a script to create random subsamples from the WCVP. It is using the
  # taxonomic and geographic information available through the data 
  # to relate species richness in the subsample to global patterns.  

# Libraries --------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(purrr)
library(progress)
library(stringr)
library(skimr)
library(sf)
library(ggplot2)
library(ggpubr)
# Defining functions -----------------------------------------------------------

# Import data ------------------------------------------------------------------
# WCVP data 
wcvp_raw <- fread("data/wcvp/wcvp_names.txt", header = T) 
dist_raw <- fread("data/wcvp/wcvp_distribution.txt", header = T) 
# TDWG regions data 
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

# manipulate data --------------------------------------------------------------

# filtering accepted species 
wcvp_accepted <- wcvp_raw %>% 
  filter(taxon_status %in% c("Accepted")) %>% 
  filter(taxon_rank %in% c("Species")) %>% 
  filter(species_hybrid == "") %>% 
  filter(genus_hybrid == "")

# filtering their native distribution
dist_native <- dist_raw %>% 
  filter(introduced == 0) %>% 
  filter(!area == "") %>% 
  filter(plant_name_id %in% wcvp_accepted$plant_name_id)

write.table(dist_native, "data/dist_native.txt")

# summarising the patterns in 3 columns 
dist_patterns <- dist_native %>%
  group_by(plant_name_id) %>%
  summarise(areas = paste(unique(area), collapse = ','),
            continent = paste(unique(continent), collapse = ','),
            region = paste(unique(region), collapse = ','))

# creating a combined dataset 
plants_full <- inner_join(wcvp_accepted, dist_patterns, by = "plant_name_id")
write.table(plants_full, "data/wcvp_accepted_merged.txt")


# extracting names from shapefile  
tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  select(LEVEL1_NAM, LEVEL1_COD)

# making it a character for later
tdwg_1_names$LEVEL1_COD <- as.character(tdwg_1$LEVEL1_COD)

# baseline richness ------------------------------------------------------------

# calculating overall richness patterns at continent level  
richness_patterns_con <- dist_native %>% 
  group_by(continent_code_l1) %>% 
  summarise(richness = n_distinct(plant_name_id))%>% 
  rename(LEVEL1_COD = continent_code_l1) %>% 
  mutate(ID = "con")

# calculating overall richness patterns at region level  
richness_patterns_reg <- dist_native %>% 
  group_by(region_code_l2) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL2_COD = region_code_l2) %>% 
  mutate(ID = "reg")

# calculating overall richness patterns at bontical recording unit level 
richness_patterns_bru <- dist_native %>% 
  group_by(area_code_l3) %>% 
  summarise(richness = n_distinct(plant_name_id)) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")

# renaming columns 
richness_patterns1 <- richness_patterns_con %>% 
  rename( LEVEL_COD = LEVEL1_COD)
richness_patterns2 <- richness_patterns_reg %>% 
  rename( LEVEL_COD = LEVEL2_COD)
richness_patterns3 <- richness_patterns_bru %>% 
  rename( LEVEL_COD = LEVEL3_COD)

#combining the datasets, keeping all rows 
richness_patterns <- rbind(richness_patterns1,richness_patterns2, 
                           richness_patterns3)

write.csv(richness_patterns, "data/richness_patterns.txt")
# Subsampling ------------------------------------------------------------------

# creating objects for subsample process 

# names 
plantlist_names <- plants_full %>%  
  select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name)

# distribution
plantlist_dist <- dist_native %>% 
  select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3)

# Region names 
tdwg_codes <- as.data.frame(tdwg_3) %>% 
  select(LEVEL3_COD,LEVEL2_COD,LEVEL1_COD)

# Overall patterns at botanical recording unit scale  
rich_overall_bru <- richness_patterns %>% 
  filter(ID =="bru") %>% 
  left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))

# Subsampling process
start.time <- Sys.time() # time keeping 

# empty objects to fill 
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
      # stepwise increase of sample size (cuts down running time)
      if (nrow(cumulative_namelist) < 10) {
        species_sample <- sample_n(plantlist_names_left, 1)
        }
        if (nrow(cumulative_namelist) > 10 & nrow(cumulative_namelist) < 100) {
          species_sample <- sample_n(plantlist_names_left, 10)
          }
          if (nrow(cumulative_namelist) > 100 & nrow(cumulative_namelist) < 1000){
            species_sample <- sample_n(plantlist_names_left, 100)
           } 
           if (nrow(cumulative_namelist) > 1000 & nrow(plantlist_names_left) > 1000){
             species_sample <- sample_n(plantlist_names_left, 1000)
              } 
              if (nrow(plantlist_names_left) < 1000){
                species_sample <- plantlist_names_left
              }
                if (nrow(plantlist_names_left) == 0){
                  break
                }
      # create cumulative dataframe with all names that have been sample
      cumulative_namelist <- rbind(cumulative_namelist, species_sample)
      species <- cumulative_namelist
      
      # extracting distribution of sample 
      dist <- dist_native %>% 
        filter(plant_name_id %in% species$plant_name_id)
      
      # richness patterns across botanical recording units 
      sample_rich_bru <- dist %>% 
        group_by(area_code_l3) %>% 
        summarise(richness_sample = n_distinct(plant_name_id)) %>% 
        rename(LEVEL_COD = area_code_l3)
      
      # combing with datasets with overall patterns 
      rich_rel <- rich_overall_bru %>% 
        left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
        replace(is.na(.), 0) %>% # replacing NA with 0 (no correlation)
        mutate(sp = nrow(cumulative_namelist))
      
      rich_cumulative <- rbind(rich_cumulative, rich_rel)
      
      # measuring correlation between richness patterns 
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
        
        # combining it all
        spear <- as.data.frame(cor.sp) %>% 
          mutate(id = rownames(.))
        rich_rel_con <- as.data.frame(cor.pe) %>% 
        mutate(id = rownames(.),
               sp = nrow(cumulative_namelist)) %>% 
        left_join(spear, by = "id")
        
        # cumulative patterns
        rich_rel_cumulative <- rbind(rich_rel_cumulative,rich_rel_bru, 
                                     rich_rel_con)
    if (nrow(cumulative_namelist) == nrow(plantlist_names)){
        break # stop when all names have been sampled 
      }
    if (nrow(cumulative_namelist) %in% seq(1,nrow(plantlist_names),10000)){
    print(paste0("There are ", nrow(cumulative_namelist) ,
                 " species in the subsample")) # print subsample size 
      }
    }
  # extracting results to objects saved outside
  samplelist[[x]] <- rich_rel_cumulative 
  samples[[x]] <- rich_cumulative
  cumulative_namelist <- data.frame()
  rich_rel_cumulative <- data.frame()
  rich_cumulative <- data.frame()
  print(paste0("This is iteration ", x))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}

# Formatting output ------------------------------------------------------------

data_out <- as.data.frame(do.call(rbind, samplelist))
data_out2 <- as.data.frame(do.call(rbind, samples))
write.table(data_out, "full_samples_rel_steps_increasing_niter100.txt")
write.table(data_out2, "full_samples_steps_increasing_niter100.txt")

# old loop 
repeat{
  
  # randomly sampling the data for species ids 
  plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist$plant_name_id)
  species_sample <- sample_n(plantlist_names_left, 1) 
  cumulative_namelist <- rbind(cumulative_namelist, species_sample)
  species <- cumulative_namelist
  # caluclating sample species richness
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # make richness patterns with brus and group by continent/region instead of ID to look at within patterns
  sample_rich_bru <- dist %>% 
    group_by(area_code_l3) %>% 
    summarise(richness_sample = n_distinct(plant_name_id)) %>% 
    rename(LEVEL_COD = area_code_l3) %>% 
    left_join(tdwg_codes, by =c("LEVEL_COD" = "LEVEL3_COD"))
  # measuring correlation
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = "LEVEL_COD") %>% 
    replace(is.na(.), 0)
  
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor(richness_sample, richness, method="spearman"),
              cor.pe = cor(richness_sample, richness, method="pearson")) %>% 
    mutate(sp = nrow(cumulative_namelist))
  
  rich_rel_cumulative <- rbind(rich_rel_cumulative,rich_rel_bru)
  # data entry
  # can try an make it into long dataframe with scale (con,reg,biu) as group 
  # variable and then calculate the means 
  
  if (nrow(cumulative_namelist) == nrow(plantlist_names)){
    break
  }
  print(paste0("There are ", nrow(cumulative_namelist) ," species in the subsample"))
}

# useful tips ------------------------------------------------------------------
# go to Edit/Folding/Collapse all to collapse all sections
# save a pdf and a png of your files
# objects names -> x_y 
# function names -> x.y
# files -> x_y_z_date.R
# for inline comments -> <space><space>#<space>comment
# library("formatR") -> function tidy_source() & tidy_dir() to clean old scripts
# to replace -> gsub(".", "_", names(dataframe), fixed = TRUE)
# change case -> tolower(names(dataframe))
