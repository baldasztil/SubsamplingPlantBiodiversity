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

options(ragg.max_dim = 10000000000)
# Defining functions -----------------------------------------------------------
theme_default.lb <- function () {
  theme(
    # add border 1)
    #panel.border = element_rect(colour = "black", fill = NA),
    # color background 2)
    panel.background = element_rect(fill = NA),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "black", linetype = 3, size = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "grey50", linetype = 3, size = 0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", face = "plain", family = "sans", size = 9),
    axis.title = element_text(colour = "black", face = "bold", family = "sans", size = 14),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black", linetype = 1, size = 0.5),
    # legend at the bottom 6)
    legend.position = "right", 
    legend.title = element_text(colour = "black", face = "plain", family = "sans", size = 14), 
    legend.text = element_text(colour = "black", face = "plain", family = "sans", size = 12), 
    
    strip.background = element_rect(colour = "black", fill = alpha("darkorange",0.5), linetype = 1, linewidth = 0.6),
    strip.text = element_text(colour = "black",face = "plain", family = "sans", size = 14),
    panel.spacing.y = unit(1, "lines"),
    panel.ontop = F)
}
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
  
  dist <- plantlist_dist %>% 
    filter(plant_name_id %in% species$plant_name_id)
  
  # phy_div patterns across tdwg3 areas
  area_codes <- unique(dist$area)
  sample_rich_bru <- data.frame(matrix(nrow = length(area_codes), ncol = 2))
  names(sample_rich_bru) <- c("area_code_l3", "phy_div_sample")
  
  
  # could also have phylogenies per country
  drop_tree_sample <-drop.tip(output_tree,no_tip_needed$label)
  
  for (i in 1: length(area_codes)) {
    temp1 <- dist %>% 
      filter(area == area_codes[i])
    
    bru_tips <- temp1 %>% 
      filter(label %in% drop_tree_sample$tip.label)
    
    sample_rich_bru$area_code_l3[i] <- temp1$area_code_l3[1]
    sample_rich_bru$phy_div_sample[i] <- pd.calc(megamatrix, tip.subset = bru_tips, method = "TBL", root.edge=FALSE)
    
    #print(paste0("Currently located in ",area_codes[i], " (",i,"/",length(area_codes),")", "|", time.laps))
  }
  
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
output_tree <- read.tree("data/phylomaker/output_tree.tre")
phy_phy_div_patterns <- fread("data/phylomaker/phylo_div_named.txt")
tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")


megamatrix <- clade.matrix(output_tree)


plants_full_accepted <- fread("data/wcvp_accepted_merged.txt") %>% 
  mutate(label = paste0(genus,"_",species))
plant_labels <- dplyr::select(plants_full_accepted, plant_name_id, label)

dist_native <- fread("data/dist_native.txt") %>% 
  left_join(plant_labels, by = "plant_name_id")
# Import data ------------------------------------------------------------------

plantlist_dist <- dist_native %>% 
 dplyr::select(plant_name_id,continent_code_l1,region_code_l2,area_code_l3, label, area)

plantlist_names <- plants_full_accepted %>% 
 dplyr::select(plant_name_id, taxon_rank, family, lifeform_description,taxon_name, label)

rich_overall_bru <- phy_phy_div_patterns 

#### Analysis 
for (i in 1:100)  {
  samples <- list()
  
  print(paste0("This is iteration ", i))
  samples[[1]] <- mclapply(seq(0,nrow(plantlist_names),1000), mc.cores = 10, subsampling.plants)
  xx <- do.call(bind_rows, samples[[1]])
  write.table(xx, paste0("data/phylogeny/Samples_iteration_phylo_",i,".txt"))
}


