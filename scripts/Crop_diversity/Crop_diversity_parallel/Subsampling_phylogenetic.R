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
library(ape)
library(U.PhyloMaker)
library(phytools)
library(caper)
library(dplyr)

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

# Import data ------------------------------------------------------------------
megatree <- read.tree("data/phylomaker/plant_megatree.tre")
gen_list <- read.csv("data/phylomaker/plant_genus_list.csv")

output_tree <- read.tree("data/phylomaker/output_tree.tre")
out_splist <- fread("data/phylomaker/output_splist.txt")


# Import data ------------------------------------------------------------------



tdwg_3 <- st_read(dsn = "data/wgsrpd-master/level3")
tdwg_2 <- st_read(dsn = "data/wgsrpd-master/level2")
tdwg_1 <- st_read(dsn = "data/wgsrpd-master/level1")
richness_patterns <- fread("data/richness_patterns.txt")
plants_full_accepted <- fread("data/wcvp_accepted_merged.txt") %>% 
  mutate(label = paste0(genus,"_",species))


missing_sp <- out_splist %>% 
  filter(output.note == "no insertion" ) 


missing_plants <- plants_full_accepted %>% 
  filter(taxon_name %in% missing_sp$species)

plant_labels <- dplyr::select(plants_full_accepted, plant_name_id, label)


dist_native <- fread("data/dist_native.txt") %>% 
  left_join(plant_labels, by = "plant_name_id")
area_codes <- unique(dist_native$area_code_l3)


tdwg_3_names <- as.data.frame(tdwg_3) %>% 
  left_join(tdwg_1, by = "LEVEL1_COD") %>% 
  dplyr::select(LEVEL3_NAM, LEVEL3_COD, LEVEL1_COD, LEVEL1_NAM)

tdwg_1_names <- as.data.frame(tdwg_1) %>% 
  dplyr::select(LEVEL1_NAM, LEVEL1_COD)

species_per_bru <- dist_native %>% 
  group_by(area_code_l3) %>% 
  reframe(plant_name_id = unique(plant_name_id)) %>% 
  left_join(tdwg_3_names, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  rename(LEVEL3_COD = area_code_l3) %>% 
  mutate(ID = "bru")


plants_formated_sample <- plants_full_accepted

temp <- dist_native %>% 
  filter(plant_name_id %in% plants_formated_sample$plant_name_id) 

bru_tips <-intersect(temp$label,
                     output_tree$tip.label)

no_tip_needed <-plants_full_accepted %>% 
  filter(!label %in% bru_tips)

drop_tree_sample <-drop.tip(output_tree,no_tip_needed$label)

area_codes <- unique(temp$area)
bru_phy_rich <- data.frame(matrix(nrow = length(area_codes), ncol = 2))
names(bru_phy_rich) <- c("area_code_l3", "phy_div")
start.time <- Sys.time() # time keeping 

for (i in 1: length(area_codes)) {
  temp1 <- temp %>% 
    filter(area == area_codes[i])
  
  bru_tips <- temp1 %>% 
    filter(label %in% drop_tree_sample$tip.label)

  
  no_tip_needed <- plants_formated_sample %>% 
    filter(!label %in% bru_tips$label)
  
  drop_tree_temp <-drop.tip(drop_tree_sample,no_tip_needed$label)
  bru_phy_rich$area_code_l3[i] <- temp1$area_code_l3[1]
  bru_phy_rich$phy_div[i] <- pd.calc(drop_tree_temp, tip.subset = NULL, method = "TBL", root.edge=FALSE)
 
  end.time <- Sys.time() # time keeping 
  time.laps <- end.time - start.time
  print(paste0("Currently located in ",area_codes[i], " (",i,"/",length(area_codes),")", "|", time.laps))
}
     
write.table(bru_phy_rich, "data/phylogeny/bru_phy_rich.txt")  
