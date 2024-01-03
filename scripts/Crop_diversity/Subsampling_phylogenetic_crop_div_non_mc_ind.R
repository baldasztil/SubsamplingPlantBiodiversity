# Subsampling Global Plant Biodiversity 13/02/2022
# data available @ ??

# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 13/02/2022

# Intro ------------------------------------------------------------------------
# This is a script to create random subsamples from the WCVP. It is using the
# taxonomic and geographic information available through the data 
# to relate species phy_div in the subsample to global patterns.  

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(ape)
library(phytools)
library(caper)
library(dplyr)
#library(picante)
#library(devtools)
#library(treeman)

# Defining functions -----------------------------------------------------------
# put in species number (spec_n) to generate a subsample
subsampling.plants <- function(spec_n) {
  
  # names in sample
  cumulative_namelist <- c()
  
  # accumulative subsampling 
  plantlist_names_left <- plantlist_names %>% 
    filter(!plant_name_id %in%  cumulative_namelist)
  if (nrow(plantlist_names_left) > spec_n) {
    species_sample <- sample_n(plantlist_names_left, spec_n)
    # to include last sample
  } else {
    species_sample <- plantlist_names_left
  }
  cumulative_namelist <- c(cumulative_namelist, species_sample$plant_name_id)
  
  # species in sample
  species <-plantlist_names %>% 
    filter(plant_name_id %in%  cumulative_namelist)
  
  # native distances of species in sample
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # codes for TDWG3 regions 
    area_codes <- unique(dist$area) 
  
  # generating a dataset to fill PD values in 
  sample_rich_bru <- data.frame(matrix(nrow = length(area_codes), ncol = 2)) 
  names(sample_rich_bru) <- c("area_code_l3", "phy_div_sample")
  
  # tips in sample
  bru_tips <-intersect(dist$label,
                       output_tree$tip.label)
  # tips to discard 
  no_tip_needed <-plants_full_accepted %>% 
    filter(!label %in% bru_tips)
  
  # trimmed tree 
  drop_tree_sample <-drop.tip(output_tree,no_tip_needed$label)
  
  #subsampling tree 
  for (i in 1: length(area_codes)) {
    
    # species in a certain TDWG3 region 
    temp1 <- dist %>% 
      filter(area == area_codes[i])
    
    # tips within the TDWG3 sample
    bru_tips <- temp1 %>% 
      filter(label %in% drop_tree_sample$tip.label)
    
    # tips to discard 
    no_tip_needed <- species %>% 
      filter(!label %in% bru_tips$label)
    
    # trimmed tree to calculate PD from 
    drop_tree_temp <-drop.tip(drop_tree_sample,no_tip_needed$label)
    
    # TDWG3 region code 
    sample_rich_bru$area_code_l3[i] <- temp1$area_code_l3[1]
    
    # calculate phylogenetic diversity for TDWG3 region
    sample_rich_bru$phy_div_sample[i] <- pd.calc(drop_tree_temp, tip.subset = NULL, method = "TBL", root.edge=FALSE)
    
    #print(paste0("Currently located in ",area_codes[i], " (",i,"/",length(area_codes),")", "|", time.laps))
  }
  
  # comibining sample with overall patterns 
  rich_rel <- rich_overall_bru %>% 
    left_join(sample_rich_bru, by = c("area_code_l3")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(sp = length(cumulative_namelist))
  
  # measuring correlation
  
  #correlation to overall phy_div 
  rich_rel_bru <- rich_rel %>% 
    summarise(cor.sp = cor.test(phy_div_sample, phy_div,
                                method="spearman", exact =F)[[4]],
              cor.pe = cor.test(phy_div_sample, phy_div,
                                method="pearson", exact =F)[[4]]) %>% 
    mutate(id = "overall",
           sp = length(cumulative_namelist))
  
  # per continent correlation
  cor.sp <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$phy_div_sample, .$phy_div, 
                      method="spearman", data = ., exact =F)[[4]])
  cor.pe <- rich_rel %>% 
    split(.$LEVEL1_COD) %>% 
    map_dbl(~cor.test(.$phy_div_sample, .$phy_div, 
                      method="pearson", data = ., exact =F)[[4]])  
  
  # combining the results
  spear <- as.data.frame(cor.sp) %>% 
    mutate(id = rownames(.))
  rich_rel_con <- as.data.frame(cor.pe) %>% 
    mutate(id = rownames(.),
           sp = length(cumulative_namelist)) %>% 
    left_join(spear, by = "id")
  
  
  if (length(cumulative_namelist) %in% seq(1,nrow(plantlist_names),1000)){
    print(paste0("There are ", length(cumulative_namelist) ," species in the subsample")) 
  }
  
  #extracting cumulative pattern
  
  return(bind_rows(rich_rel_bru, rich_rel_con))
}

# Import data ------------------------------------------------------------------

#phylo tree
output_tree <- read.tree("data/phylomaker/output_tree.tre")

# phylo diversity overall patterns 
phy_phy_div_patterns <- fread("data/phylomaker/phylo_div_named.txt")

# tdwg regions 
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")

# WCVP full species list 
plants_full_accepted <- fread("data/wcvp_accepted_merged.txt") %>% 
  mutate(label = paste0(genus,"_",species))
plant_labels <- dplyr::select(plants_full_accepted, plant_name_id, label)


# native ranges of species 
dist_native <- fread("data/dist_native.txt") %>% 
  left_join(plant_labels, by = "plant_name_id")
# Format data ------------------------------------------------------------------

plantlist_dist <- dist_native %>% 
 dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3, label, area)

plantlist_names <- plants_full_accepted %>% 
 dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name, label)

rich_overall_bru <- phy_phy_div_patterns 

#### Analysis 
for (i in 1:1)  {
  samples <- list() # sample patterns 
  print(paste0("This is iteration ", i)) # iteration number 
  samples[[1]] <- lapply(seq(1,nrow(plantlist_names),1000), subsampling.plants) # apply function to accumulate samples in steps of 1000, 
  #this is normally parallelised on Linux with the mclapply function but that does not work on Windows 
  xx <- do.call(bind_rows, samples[[1]]) # tranform list output to dataframe
  write.table(xx, paste0("data/phylogeny/Samples_iteration_phylo_",i,".txt")) # safe results as txt file
}


